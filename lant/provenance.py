"""What produced a result: code version, environment, third-party programs; a lock so that one output directory has
one writer at a time.

  <output>/.lant/provenance/<id>.json   written once per distinct (code, environment); every entry of run_log.jsonl
                                          names its id.  Code: git commit, whether the working tree differs from it
                                          and a fingerprint (sha1) of that difference and of untracked source files;
                                          Python and package versions; platform; sha1 of the third-party programs.
  <output>/.lant/lock                   the process that is writing to the output directory (PID, host, start).
"""
import datetime as _dt
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGES = ("numpy", "scipy", "h5py", "obspy", "matplotlib", "yaml", "openpyxl")


def _sha1_file(path, chunk=1 << 20):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def _git(*args):
    try:
        r = subprocess.run(["git", "-C", ROOT] + list(args), capture_output=True, text=True, timeout=30)
        return r.stdout if r.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def code_version():
    head = _git("rev-parse", "HEAD")
    if head is None:
        from . import __version__
        return dict(git=False, version=__version__)
    diff = _git("diff", "HEAD") or ""
    untracked = [f for f in (_git("ls-files", "--others", "--exclude-standard") or "").split("\n")
                 if f.endswith((".py", ".f90", ".yaml", ".sh", ".patch"))]
    return dict(git=True, commit=head.strip(), modified=bool(diff), diff_sha1=hashlib.sha1(diff.encode()).hexdigest() if diff else None,
                untracked_sources={f: _sha1_file(os.path.join(ROOT, f)) for f in sorted(untracked)
                                   if os.path.isfile(os.path.join(ROOT, f))})


def environment():
    from .tools import external_path
    pk = {}
    for m in PACKAGES:
        try:
            pk[m] = getattr(__import__(m), "__version__", "?")
        except ImportError:
            pk[m] = None
    ext = {}
    for b in ("DSurfTomo", "pDSurfTomo", "fwd1d"):
        p = external_path("bin", b)
        if os.path.isfile(p):
            ext[b] = dict(path=p, sha1=_sha1_file(p))
    return dict(python=sys.version.split()[0], packages=pk, platform=platform.platform(), machine=platform.machine(),
                external=ext)


def record(rec_dir):
    """Write .lant/provenance/<id>.json (once) and return the id."""
    info = dict(code=code_version(), environment=environment())
    pid = hashlib.sha1(json.dumps(info, sort_keys=True).encode()).hexdigest()[:12]
    d = os.path.join(rec_dir, "provenance")
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "%s.json" % pid)
    if not os.path.exists(p):
        info["first_used"] = _dt.datetime.now().isoformat(timespec="seconds")
        with open(p + ".tmp", "w") as fh:
            json.dump(info, fh, indent=1)
        os.replace(p + ".tmp", p)
    return pid


class Lock:
    """.lant/lock: one writing process per output directory, held with an exclusive flock (released by the
    system when the process ends, also when it is killed).  The file says which process holds it.  A nested step
    of the same process passes."""
    _held = {}

    def __init__(self, rec_dir):
        self.path = os.path.join(rec_dir, "lock")
        self.fh = None

    def __enter__(self):
        import fcntl
        if self.path in Lock._held:                  # nested step of the same process
            Lock._held[self.path][1] += 1
            return self
        fh = open(self.path, "a+")
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            fh.seek(0)
            other = fh.read().strip()
            fh.close()
            raise RuntimeError("another lant process is writing to this output directory (%s); wait for it to finish"
                               % (other or "unknown process"))
        fh.seek(0)
        fh.truncate()
        fh.write(json.dumps(dict(pid=os.getpid(), host=socket.gethostname(),
                                 start=_dt.datetime.now().isoformat(timespec="seconds"))))
        fh.flush()
        self.fh = fh
        Lock._held[self.path] = [fh, 1]
        return self

    def __exit__(self, *a):
        import fcntl
        ent = Lock._held.get(self.path)
        if ent is None:
            return False
        ent[1] -= 1
        if ent[1] == 0:
            fh = ent[0]
            fh.seek(0)
            fh.truncate()
            fcntl.flock(fh, fcntl.LOCK_UN)
            fh.close()
            del Lock._held[self.path]
        return False

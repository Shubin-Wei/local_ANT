"""Confirmation points: the chain stops only where a scientific judgement is the user's to make.

After a step, lant.confirm.judgements lists what the program cannot decide (lant.manifest keeps the list and a
fingerprint, the "basis", of the results it is about).  When the list is empty the chain goes on.  Otherwise
`lant confirm STEP -c CONFIG` writes the confirmation into the `confirmed:` block at the end of the configuration
file, between two marker lines (the rest of the file is not touched), with the time and the basis.  The
configuration file then holds every decision: a rerun with the same results does not stop again, and when the
results change the basis no longer matches and the step has to be confirmed again.

What is judged at each point (parameters that follow mechanically from the data are not confirmed; they are
listed with their derivation in .lant/steps/<step>.md):
  clean      station names that could not be matched one to one, coordinates differing from the record headers,
             channels whose orientation the names do not tell, stations with several sensors
  stack      the candidate station-days / station-hours, when there are any; nothing is removed unless confirmed
             with --exclude-candidates
  reference  whether the F-J picks follow the fundamental mode (a borrowed curve is not confirmed)
  inversion  the starting model and its fit, the depth nodes (deepest node) and the velocity bounds; the accepted
             depth nodes and bounds are written into the configuration

Precedence of a value: a number written by the user in the configuration > a confirmed value (while its
confirmation holds) > the rule's value.  "auto" written in the configuration is the same as leaving the key out.
"""
import datetime as _dt
import hashlib
import json
import os

import numpy as np
import yaml

from . import auto, paths
from . import config as C
CONFIRM_STEPS = C.CONFIRM_STEPS          # used by callers as confirm.CONFIRM_STEPS

BEGIN = "# ==== confirmed: written by `lant confirm` (values accepted at the confirmation points; edit the values if needed, keep the marker lines) ===="
END = "# ==== end of confirmed ===="
VALUE_KEYS = {"stack": ["stack.exclude_candidates"],             # the values a confirmation carries
              "inversion": ["inversion.depth_sequence_km", "inversion.vmin", "inversion.vmax"]}


class NeedsConfirmation(RuntimeError):
    pass


def _hash(*items):
    h = hashlib.sha1()
    for it in items:
        if isinstance(it, str) and os.path.isfile(it):
            with open(it, "rb") as fh:
                h.update(fh.read())
        else:
            h.update(json.dumps(it, sort_keys=True, default=str).encode())
    return h.hexdigest()[:16]


def _rule_values(cfg, keys):
    dec = auto.load(cfg).get("parameters", {})
    return {k: dec[k]["value"] for k in keys if k in dec}


def proposals(cfg, step):
    """{dotted key: value} that `lant confirm step` writes into the configuration: the rule's values of the keys
    the user did not give."""
    return {k: v for k, v in _rule_values(cfg, VALUE_KEYS.get(step, []) if step != "stack" else []).items()
            if auto.user_value(cfg, k.split(".")) is None}


def judgements(cfg, step, summary):
    """The scientific judgements `step` leaves to the user (texts); empty when there is nothing to judge."""
    s = summary or {}
    out = []
    if step == "clean":
        out += [n for n in s.get("channel_notes", []) if "confirm" in n or "not used" in n]
        names = s.get("names") or {}
        for key, what in (("ambiguous", "station names matching more than one table entry"),
                          ("table_without_data", "stations of the table without data"),
                          ("data_without_table", "records without a station in the table"),
                          ("coordinate_mismatch_km", "coordinates differing from the record headers (km)")):
            v = names.get(key)
            if v:
                out.append("%s: %s" % (what, _short(v)))
        at = s.get("array_type") or {}
        if not at.get("given") and at.get("clear") is False:
            out.append("array type taken as %s: %s.  If it is the other type, write project.array_type in the "
                       "configuration" % (at.get("type"), at.get("text")))
        r = s.get("response") or {}
        if r and not r.get("removed") and not (len(r.get("models") or {}) == 1 and r.get("named") == r.get("stations")):
            out.append("instrument responses are not removed (no data.response): are all stations the same "
                       "instrument model?  If not, give the response files in data.response")
        if s.get("several_sensors"):
            out.append("stations with several sensors or rates (one kept): %s" % _short(sorted(s["several_sensors"])))
    elif step == "stack":
        c = s.get("candidates") or {}
        n = c.get("station_blocks", 0)
        if n > c.get("expected_by_chance", 0):       # as many as chance alone gives: listed in the report, no stop
            out.append("%d station sub-stacks differ from the rest (03_stack/candidates.txt, station_qc/station_qc.png)"
                       ": remove them (`lant confirm stack --exclude-candidates`) or keep them (`lant confirm stack`)?"
                       % n)
    elif step == "reference":
        if s.get("source") == "fj":
            out.append("do the F-J picks follow the fundamental-mode ridge (04_reference/fj.png)?  Only the fundamental "
                       "Rayleigh mode is measured; the strongest ridge can be a higher mode")
            if s.get("zones"):
                out.append("zone references: %d of %d zones have their own curve (04_reference/zones.png, zones.txt): "
                           "do the zone curves differ as the geology suggests, and none jumps to another ridge?"
                           % (s["zones"]["with_own_curve"], s["zones"]["zones"]))
            if cfg["project"].get("array_type") == "line" and \
                    list(cfg["data"]["velocity_range_kms"]) == list(C.DEFAULTS["data"]["velocity_range_kms"]):
                out.append("data.velocity_range_kms is the default %s on a line: the picks may then cover a narrow "
                           "band only (fj.txt), and beyond it faster arrivals can be accepted in the disp step.  If "
                           "the velocities of the area are known, write them and run again"
                           % list(cfg["data"]["velocity_range_kms"]))
            out += list(s.get("warnings") or [])
            if s.get("held_flat"):
                out.append(s["held_flat"])
    elif step == "inversion":
        st = s.get("start") or {}
        seq = cfg["inversion"]["depth_sequence_km"]
        out.append("starting model and its fit to the average curve (06_inversion/vs1d_fit.txt, figures/4_start_model): "
                   "acceptable?")
        if isinstance(seq, (list, tuple)) and seq:
            out.append("depth nodes down to %.3g km (%d nodes): deep enough, and not deeper than the data constrain?"
                       % (float(seq[-1]), len(seq)))
        out.append("velocity bounds %s-%s km/s: do they bracket the Vs expected in the area?"
                   % (st.get("vmin", cfg["inversion"]["vmin"]), st.get("vmax", cfg["inversion"]["vmax"])))
    return out


def _short(v, n=160):
    t = json.dumps(v, ensure_ascii=False, default=str) if not isinstance(v, str) else v
    return t if len(t) <= n else t[:n - 3] + "..."


def basis(cfg, step, items):
    """Fingerprint of the results the judgements `items` are about (None when there is nothing to judge)."""
    if not items:
        return None
    P = lambda *p: paths.P(cfg, *p)
    if step == "clean":
        return _hash(items)
    if step == "stack":
        return _hash(P(paths.STACK, "candidates.txt"))
    if step == "reference":
        z = np.load(P(paths.REFERENCE, "reference.npz"))
        zp = P(paths.REFERENCE, "zones.npz")
        zones = np.round(np.load(zp)["c_ref"], 4).tolist() if os.path.exists(zp) else None
        return _hash(np.round(z["freqs"], 6).tolist(), np.round(np.where(z["ref_valid"], z["c_ref"], 0), 4).tolist(), zones)
    if step == "inversion":
        return _hash(P(paths.INV, "average_curve.txt"), P(paths.INV, "vs1d.txt"),
                     _rule_values(cfg, VALUE_KEYS["inversion"]))
    raise ValueError(step)


def valid_values(cfg):
    """{dotted key: value} of the confirmations that still hold (their basis is that of the step's present
    results)."""
    from . import manifest
    steps = manifest.load(cfg)["steps"]
    out = {}
    for step, keys in VALUE_KEYS.items():
        c = (cfg.get("_confirmed") or {}).get(step) or {}
        r = steps.get(step)
        if r and r.get("basis") and c.get("_basis") == r["basis"]:
            out.update({k: c[k] for k in keys if k in c})
    return out


def _read_block(text):
    """(text without the confirmed block, the block as a dict)."""
    if BEGIN not in text:
        if "\nconfirmed:" in "\n" + text:
            raise C.ConfigError("the configuration has a confirmed: block without the marker lines of `lant confirm`; "
                                "move its values into the configuration or delete it")
        return text, {}
    head, rest = text.split(BEGIN, 1)
    if END not in rest:
        raise C.ConfigError("the confirmed block of the configuration has no end marker (%s)" % END)
    block, tail = rest.split(END, 1)
    data = yaml.safe_load(block) or {}
    return head.rstrip("\n") + "\n" + tail.lstrip("\n"), data.get("confirmed", {})


def _write_block(path, conf):
    text = open(path).read()
    body, _ = _read_block(text)
    dumped = yaml.safe_dump({"confirmed": conf}, sort_keys=False, allow_unicode=True, default_flow_style=None, width=110)
    new = body.rstrip("\n") + "\n\n" + BEGIN + "\n" + dumped + END + "\n"
    tmp = path + ".tmp"
    open(tmp, "w").write(new)
    os.replace(tmp, path)


def confirm(cfg, step, exclude_candidates=False):
    """Write the confirmation of `step` into the configuration file -> (entry, steps to run again)."""
    from . import manifest
    if step not in C.CONFIRM_STEPS:
        raise ValueError("steps to confirm: %s" % ", ".join(C.CONFIRM_STEPS))
    r = manifest.load(cfg)["steps"].get(step)
    if r is None or r["state"] != "done":
        raise RuntimeError("step %s has not finished: `lant run all -c %s` first" % (step, cfg["_config_file"]))
    if not r.get("to_confirm") and not (step == "stack" and exclude_candidates):
        raise RuntimeError("step %s leaves nothing to confirm: go on with `lant run all -c %s`"
                           % (step, cfg["_config_file"]))
    path = cfg["_config_file"]
    _, conf = _read_block(open(path).read())
    entry = {"_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"), "_basis": r.get("basis"),
             "_judged": list(r.get("to_confirm") or [])}
    entry.update(proposals(cfg, step))
    rerun = []
    if step == "stack":
        if exclude_candidates:
            entry["stack.exclude_candidates"] = True
            rerun.append("stack")
        else:
            entry["_decision"] = "candidates kept (nothing removed)"
    conf[step] = entry
    _write_block(path, conf)
    return entry, rerun

"""Main-process progress on stderr; unknown work uses a heartbeat, never a fake ETA."""

import os
import sys
import threading
import time

_current = None


class Progress:
    def __init__(self, label, stream=None, interval=None):
        self.label = label
        self.stream = stream or sys.stderr
        self.tty = self.stream.isatty()
        self.interval = interval if interval is not None else (1 if self.tty else 30)
        self.disabled = os.environ.get("LANT_PROGRESS", "auto").lower() in ("off", "0", "false")
        self.done, self.total, self.unit = 0, None, "tasks"
        self.detail = ""
        self.started = time.monotonic()
        self.stop = threading.Event()
        self.guard = threading.Lock()

    def render(self):
        if self.disabled:
            return
        with self.guard:
            elapsed = time.monotonic() - self.started
            count = "%d/%d %s" % (self.done, self.total, self.unit) if self.total is not None else "running"
            eta = ""
            if self.total and self.done and elapsed >= 2:
                eta = ", ETA ~%.0fs" % (elapsed * max(0, self.total - self.done) / self.done)
            bar = ""
            if self.tty and self.total:
                n = min(20, int(20 * self.done / self.total))
                bar = "[" + "=" * n + " " * (20 - n) + "] "
            line = "[lant] %s %s%s, elapsed %.0fs%s" % (self.label, bar, count, elapsed, eta)
            if self.detail:
                line += "; " + self.detail
            self.stream.write(("\r\033[K" if self.tty else "") + line + ("" if self.tty else "\n"))
            self.stream.flush()

    def _heartbeat(self):
        while not self.stop.wait(self.interval):
            self.render()

    def __enter__(self):
        global _current
        self.previous = _current
        _current = self
        self.render()
        self.thread = threading.Thread(target=self._heartbeat, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *exc):
        global _current
        self.stop.set()
        self.thread.join()
        self.render()
        if self.tty and not self.disabled:
            self.stream.write("\n")
        _current = self.previous


def track(items, total, unit):
    """Count completed tasks without changing their order or scheduling."""
    p = _current
    if p:
        with p.guard:
            p.done, p.total, p.unit = 0, total, unit
    for item in items:
        yield item
        if p:
            with p.guard:
                p.done += 1


def detail(text):
    """Observed substage information; it does not imply a percent or reliable ETA."""
    if _current:
        with _current.guard:
            _current.detail = text

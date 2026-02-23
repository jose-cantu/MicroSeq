from __future__ import annotations

import os
import platform


def _select_qt_platform() -> None:
    """Select Qt backend before any Qt/PySide imports happen."""
    # Keeping commentary context here (moved from main_window.py) so future edits stay readable like the rest of this project.
    # Goal remains the same: pick the first visible backend that has a server and avoid Wayland maximize crashes.

    # Respect explicit Qt platform from user/environment first.
    if os.environ.get("QT_QPA_PLATFORM"):
        os.environ["MICROSEQ_QT_BACKEND_REASON"] = "explicit_qt_qpa_platform"
        return

    override = (os.environ.get("MICROSEQ_QT_BACKEND") or "").strip().lower()
    if override in {"wayland", "xcb", "offscreen"}:
        os.environ["QT_QPA_PLATFORM"] = override
        os.environ["MICROSEQ_QT_BACKEND_REASON"] = "microseq_qt_backend_override"
        return

    # Linux (including WSL) backend choice should follow active display sockets.
    # This intentionally avoids relying on XDG_SESSION_TYPE because in real user shells that value can be blank.
    wayland = bool(os.environ.get("WAYLAND_DISPLAY"))
    x11 = bool(os.environ.get("DISPLAY"))

    if platform.system() == "Linux":
        if wayland and x11:
            # Stability first: Wayland session with Xwayland available defaults to xcb.
            os.environ["QT_QPA_PLATFORM"] = "xcb"
            os.environ["MICROSEQ_QT_BACKEND_REASON"] = "linux_wayland_with_xwayland_default_xcb"
        elif wayland:
            # Native wayland path when no DISPLAY/Xwayland is available.
            os.environ["QT_QPA_PLATFORM"] = "wayland"
            os.environ["MICROSEQ_QT_BACKEND_REASON"] = "linux_wayland_no_display_default_wayland"
        elif x11:
            # Plain X11 session path.
            os.environ["QT_QPA_PLATFORM"] = "xcb"
            os.environ["MICROSEQ_QT_BACKEND_REASON"] = "linux_x11_default_xcb"
        else:
            # Headless CI / SSH / no visible display server.
            os.environ["QT_QPA_PLATFORM"] = "offscreen"
            os.environ["MICROSEQ_QT_BACKEND_REASON"] = "headless_default_offscreen"

    # macOS defaults to cocoa so no override needed here unless user opts in explicitly.


def launch() -> int:
    _select_qt_platform()
    from .main_window import launch as _launch

    return _launch()

if __name__ == "__main__":
    raise SystemExit(launch())

"""Standalone PySide6 desktop GUI for MeltScope (WinMelt-style workbench).

The Qt modules are imported lazily, so data-only submodules such as
``varmelt.gui.examples`` stay importable on machines without PySide6.
"""

__all__ = ["MainWindow", "build_app", "main"]


def __getattr__(name):
    if name in __all__:
        from . import main as _main
        return getattr(_main, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

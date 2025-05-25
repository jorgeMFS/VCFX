import subprocess
import shutil
import functools

__all__ = ["available_tools", "run_tool"]


def available_tools():
    """Return a list of VCFX tools available on the PATH."""
    result = subprocess.run(["vcfx", "--list"], capture_output=True, text=True)
    if result.returncode != 0:
        return []
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def run_tool(tool, *args, check=True, capture_output=False, text=True, **kwargs):
    """Run a VCFX tool using subprocess.

    Parameters
    ----------
    tool : str
        Name of the tool without the ``VCFX_`` prefix.
    *args : list
        Arguments passed to the tool.
    check : bool, optional
        If ``True`` (default) raise ``CalledProcessError`` on non-zero
        return code.
    capture_output : bool, optional
        If ``True`` capture stdout/stderr and return them on the returned
        ``CompletedProcess`` object.
    text : bool, optional
        If ``True`` decode output as text. Defaults to ``True``.
    **kwargs : dict
        Additional keyword arguments forwarded to ``subprocess.run``.

    Returns
    -------
    subprocess.CompletedProcess
    """
    exe = shutil.which(f"VCFX_{tool}")
    if exe is None:
        raise FileNotFoundError(f"VCFX tool '{tool}' not found in PATH")
    cmd = [exe, *map(str, args)]
    return subprocess.run(cmd, check=check, capture_output=capture_output, text=text, **kwargs)


# Lazy attribute access for tool wrappers

def __getattr__(name):
    if name in available_tools():
        return functools.partial(run_tool, name)
    raise AttributeError(f"module 'vcfx' has no attribute '{name}'")

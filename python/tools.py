import subprocess
import shutil
import functools
import csv
import os
from typing import Any, Callable, Sequence

# Cache for storing the list of available tools once discovered
_TOOL_CACHE: list[str] | None = None

__all__ = [
    "available_tools",
    "run_tool",
    "alignment_checker",
    "allele_counter",
    "variant_counter",
]


def available_tools(refresh: bool = False) -> list[str]:
    """Return the list of available VCFX command line tools.

    Parameters
    ----------
    refresh : bool, optional
        If ``True`` ignore any cached value and re-run ``vcfx --list``.
        Defaults to ``False``.

    Returns
    -------
    list[str]
        Names of tools discovered on ``PATH``.

    Raises
    ------
    FileNotFoundError
        If the ``vcfx`` executable cannot be found.
    """
    global _TOOL_CACHE
    if _TOOL_CACHE is not None and not refresh:
        return _TOOL_CACHE

    exe = shutil.which("vcfx")
    if exe is None:
        tools: set[str] = set()
        for path in os.environ.get("PATH", "").split(os.pathsep):
            if not path:
                continue
            try:
                for entry in os.listdir(path):
                    if entry.startswith("VCFX_"):
                        full = os.path.join(path, entry)
                        if os.path.isfile(full) and os.access(full, os.X_OK):
                            tools.add(entry[5:])
            except OSError:
                continue
        if not tools:
            raise FileNotFoundError("vcfx wrapper not found in PATH")
        _TOOL_CACHE = sorted(tools)
        return _TOOL_CACHE

    result = subprocess.run([exe, "--list"], capture_output=True, text=True)
    if result.returncode != 0:
        _TOOL_CACHE = []
    else:
        _TOOL_CACHE = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    return _TOOL_CACHE


def run_tool(
    tool: str,
    *args: str,
    check: bool = True,
    capture_output: bool = False,
    text: bool = True,
    **kwargs: Any,
) -> subprocess.CompletedProcess:
    """Run a VCFX tool using :func:`subprocess.run`.

    Parameters
    ----------
    tool : str
        Name of the tool without the ``VCFX_`` prefix.
    *args : str
        Command line arguments passed to the tool.
    check : bool, optional
        If ``True`` (default) raise ``CalledProcessError`` on a non-zero
        exit status.
    capture_output : bool, optional
        Capture standard output and error and attach them to the returned
        ``CompletedProcess``. Defaults to ``False``.
    text : bool, optional
        If ``True`` decode output as text. Defaults to ``True``.
    **kwargs : Any
        Additional keyword arguments forwarded to :func:`subprocess.run`.

    Returns
    -------
    subprocess.CompletedProcess
        The completed process instance for the invoked command.

    Raises
    ------
    FileNotFoundError
        If the requested tool cannot be found on ``PATH``.
    subprocess.CalledProcessError
        If ``check`` is ``True`` and the process exits with a non-zero
        status.
    """
    exe = shutil.which(f"VCFX_{tool}")
    if exe is None:
        raise FileNotFoundError(f"VCFX tool '{tool}' not found in PATH")
    cmd = [exe, *map(str, args)]
    return subprocess.run(cmd, check=check, capture_output=capture_output, text=text, **kwargs)


def alignment_checker(vcf_file: str, reference: str) -> list[dict]:
    """Run ``alignment_checker`` and parse the TSV output.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file to check.
    reference : str
        Reference FASTA file.

    Returns
    -------
    list[dict]
        Parsed rows from the discrepancy report.
    """

    result = run_tool(
        "alignment_checker",
        "--alignment-discrepancy",
        vcf_file,
        reference,
        capture_output=True,
        text=True,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def allele_counter(vcf_file: str, samples: Sequence[str] | None = None) -> list[dict]:
    """Run ``allele_counter`` and return allele counts.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF input file.
    samples : Sequence[str] | None, optional
        Optional subset of sample names to process.

    Returns
    -------
    list[dict]
        Parsed allele counts.
    """

    args: list[str] = []
    if samples:
        args.extend(["--samples", " ".join(samples)])

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "allele_counter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    reader = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
    return list(reader)


def variant_counter(vcf_file: str, strict: bool = False) -> int:
    """Count variants in a VCF using ``variant_counter``.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file.
    strict : bool, optional
        If ``True`` fail on malformed lines. Defaults to ``False``.

    Returns
    -------
    int
        Number of valid variants reported by the tool.
    """

    args: list[str] = []
    if strict:
        args.append("--strict")

    with open(vcf_file, "r", encoding="utf-8") as fh:
        inp = fh.read()

    result = run_tool(
        "variant_counter",
        *args,
        capture_output=True,
        text=True,
        input=inp,
    )

    out = result.stdout.strip()
    if ":" in out:
        out = out.split(":", 1)[1]
    return int(out.strip())


# Lazy attribute access for tool wrappers

def __getattr__(name: str) -> Callable[..., subprocess.CompletedProcess]:
    """Return a callable wrapper for a VCFX tool.

    Parameters
    ----------
    name : str
        Name of the tool as exposed on ``PATH`` without the ``VCFX_`` prefix.

    Returns
    -------
    Callable[..., subprocess.CompletedProcess]
        A function that invokes the requested tool.

    Raises
    ------
    AttributeError
        If *name* does not correspond to an available tool.
    FileNotFoundError
        If the ``vcfx`` wrapper cannot be located on ``PATH``.
    """
    tools = available_tools()
    if name in tools:
        return functools.partial(run_tool, name)
    raise AttributeError(f"module 'vcfx' has no attribute '{name}'")


"""
Lazy dispatch from Click commands to the implementing Python functions.

Every Click subcommand collects its parsed options into a flat ``params`` dict
and hands it to :func:`run` together with the dotted ``module`` path and the
``func`` name that implements it.  :func:`run`

1. translates CLI-only conveniences (``--strand``) into the boolean parameters
   the API actually takes,
2. validates CRAM/reference requirements and contig compatibility, then
3. lazily imports ``module``, looks up ``func``, filters ``params`` down to the
   function's signature, and calls it.

Keeping the heavy imports (``pysam``, the ``frag`` implementations) inside
:func:`run` means ``import finaletoolkit.cli`` and ``--help`` stay fast.
"""
from __future__ import annotations

import importlib
from inspect import getfullargspec
from sys import stderr
from typing import Any


def _translate_strand(params: dict[str, Any]) -> None:
    """Map ``--strand {both,forward,reverse}`` onto the two API booleans.

    The motif functions take ``both_strands`` and ``negative_strand``.  Exposing
    those directly produced a confusing ``Default: True`` toggle, so the CLI
    offers a single ``--strand`` choice and we expand it here.
    """
    if "strand" not in params:
        return
    strand = params.pop("strand")
    params["both_strands"] = strand == "both"
    params["negative_strand"] = strand == "reverse"


def _validate_inputs(params: dict[str, Any]) -> None:
    """Validate CRAM reference requirements and contig compatibility.

    Mirrors the pre-Click behavior: CRAM input requires a reference, and when a
    reference accompanies a BAM/CRAM the contigs must be compatible.  Heavy
    imports stay local so unrelated subcommands don't pay for them.
    """
    input_file = params.get("input_file")
    reference_file = params.get("reference_file") or params.get("refseq_file")

    if not input_file:
        return

    lowered = str(input_file).lower()

    if lowered.endswith(".cram") and not reference_file:
        stderr.write(
            "Error: CRAM files require a reference file (-r/--reference).\n"
        )
        raise SystemExit(1)

    if reference_file and (lowered.endswith(".bam") or lowered.endswith(".cram")):
        import pysam

        from finaletoolkit.io.reference import ReferenceWrapper
        from finaletoolkit.utils.validation import validate_compatible_contigs

        try:
            try:
                with pysam.AlignmentFile(
                    input_file, "r", reference_filename=reference_file
                ) as sam:
                    input_contigs = dict(zip(sam.references, sam.lengths))
            except Exception as e:
                stderr.write(f"Error opening alignment file '{input_file}': {e}\n")
                raise SystemExit(1)

            try:
                with ReferenceWrapper(reference_file) as ref:
                    ref_contigs = ref.chroms
            except Exception as e:
                stderr.write(f"Error opening reference file '{reference_file}': {e}\n")
                raise SystemExit(1)

            validate_compatible_contigs(
                ref_contigs,
                input_contigs,
                validate_sizes=True,
                throw_on_error=True,
            )
        except ImportError:
            pass
        except (ValueError, RuntimeError) as e:
            stderr.write(f"Validation Error: {e}\n")
            raise SystemExit(1)


def run(module: str, func: str, params: dict[str, Any]) -> Any:
    """Validate, then lazily import ``module.func`` and call it with ``params``.

    ``params`` is filtered to the implementing function's signature (unless the
    function accepts ``**kwargs``), so CLI-only keys are dropped automatically.
    Each option's name already matches its API parameter, so no renaming is
    needed beyond the ``--strand`` expansion done here.
    """
    params = dict(params)
    _translate_strand(params)
    _validate_inputs(params)

    mod = importlib.import_module(module)
    function = getattr(mod, func)

    spec = getfullargspec(function)
    if spec.varkw is None:
        params = {
            key: value
            for key, value in params.items()
            if key in spec.args or key in spec.kwonlyargs
        }
    return function(**params)

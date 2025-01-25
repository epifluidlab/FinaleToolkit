Using Finaletoolkit Commands
-----------------------------

*   Finaletoolkit commands are specified in ``params.yaml`` with underscores instead of hyphens.
*   Set command flags by appending ``_<flag_name>`` to the converted command name.
*   **Dependencies:**  The workflow respects command dependencies.  For example, ``adjust-wps`` requires ``wps`` output, and ``mds`` needs ``end-motifs``.
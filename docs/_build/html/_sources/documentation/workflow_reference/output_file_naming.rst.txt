Output File Naming
------------------

*   **Filtered Files:** Files are always given a ``.filtered`` extension before the file format when passed into the output directory (e.g., ``file.filtered.bed.gz``).
*   **Command-Processed Files:** Files processed by a Finaletoolkit command have the command name (with underscores) inserted before their format (e.g., ``file.frag_length_bins.bed.gz``).
*   **Multiple Commands:** Input files will be processed for each enabled Finaletoolkit command.
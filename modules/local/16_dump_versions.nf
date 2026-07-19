// Collate every process's versions.yml into one software_versions.yml.

process DUMP_VERSIONS {
    label 'process_single'
    container "quay.io/biocontainers/python@sha256:f6b44640f06e8265ebf5ce85ca12cea53af110c188d6b4acf5f59887c24abb8f"  // python 3.10

    input:
    path versions

    output:
    path "software_versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    # Merge the per-process YAML blocks. A plain `sort -u` sorts *lines* and so
    # destroys the nesting (process headers end up interleaved with tool
    # versions), which is what the first version of this process did.
    import re
    from pathlib import Path

    merged = {}
    current = None

    for line in Path("${versions}").read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped == "END_VERSIONS":
            continue
        # Process header, e.g.  "CALL_VARIANTS:CLAIR3":
        # Leading whitespace is tolerated: Nextflow only strips the *common*
        # indentation of a script block, so a process whose script contains an
        # unindented heredoc keeps its versions block indented.
        header = re.match(r'^\\s*"([^"]+)":\\s*\$', line)
        if header:
            current = header.group(1)
            merged.setdefault(current, {})
            continue
        # Tool line, e.g.      clair3: v1.0.10
        if current and ":" in stripped:
            tool, _, version = stripped.partition(":")
            merged[current][tool.strip()] = version.strip()

    with open("software_versions.yml", "w") as fh:
        for proc in sorted(merged):
            fh.write('"%s":\\n' % proc)
            for tool in sorted(merged[proc]):
                fh.write("    %s: %s\\n" % (tool, merged[proc][tool]))
    """

    stub:
    """
    touch software_versions.yml
    """
}

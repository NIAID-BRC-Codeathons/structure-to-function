"""The full report — the pipeline's canonical deliverable.

Runs last, after every module has written what it has, and reads the whole run
directory rather than any one module's output. `s2f.m6_report` renders the compact
report from `report.json` alone; this renders everything the run produced.
"""

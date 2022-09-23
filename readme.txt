This python script reads all IR Data CSVs in a folder from a given path, normalizes/baseline-corrects them, and plots them, either as superposed scatterplots or bar graphs for a peak area within each spectrum

Resources:
 - J Hofman's 2016 "IR Spectroscopic Method for Determiation of Silicone Cross-Linking"
 - Joseph Fortenbaugh's 2019 PSU Doctoral Thesis

To-do
[General Graphing]
    - Spectra:
        - better color system (auto assign from cmap, call on local .txt, cluster into Hues for trial and move between saturation within trials, etc.)
            - clever way to group the separate categories? then call on like 'blues.txt, reds.txt, etc.' for hues within that category?
    - Bar Graph
        - plot all the regions in shared window (with the scatterplot as a sub-window?)
        - use averages and graph bars with standard error/dev first (how to procedurally take arbitrary number of values to average)
            - cluster them based on appearance of suffix '1'!
        - in integration bar graph: procedurally color bars for each set to separate without label names
        - have bar labels end on tick instead of middle on tick

[General Readability/Usability/Efficiency]
    - rewrite abs conversion and peak integration as methods
    - exclude lowest WNs that always message up the y scaling
    - drag/drop GUI (input file directory, graphs, wavenumber bounds, checkbox for included plot types and export) and web implementation'
        - directory fields (1: CSV folder, 2: output folder)
        - number field: controls to average over
        - checkboxes: Graphs (1: line, 2: bar), Bands (if bar selection, select all wanted to create)
        - buttons: 1: Run, 2: Export results CSV

[Interpretation]
    - what do first two areas represent, how should they change with curing, and why do they go opposite sometimes with laser cure?
    - thinking: homogeneity of samples in ATR?
    - path length - cb change dielectric meaningfully? (what is dielectic constant for each)
    - creating some surface effects that could lead to scatter in measurement?
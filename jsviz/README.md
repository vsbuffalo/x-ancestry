# Javascript Genealogy Visualizations

Two plots in the paper are the result of the javascript and
[d3js](http://d3js.org) visualizations. For completeness I have included the
source files here, but the code is available in a separate repository,
[vsbuffalo/genealogyjs](https://github.com/vsbuffalo/genealogyjs). The code to
generate the graphics is from the simulation code, with `--seed 3`:

    python ../ancestralsim.py x --ngens 9 --seed 3 --nsims 5 --json > x.json

Then start the Python HTTP server:

    python -m SimpleHTTPServer 8082

Then, the visualizations are http://127.0.0.1:8082/shared.html (shared
segments) and http://127.0.0.1:8082/ (family arc). Writing SVGs from d3 is
tricky, so I use the New York Times terrific tool, [SVG
Crowbar](http://nytimes.github.io/svg-crowbar/) to pull down the `.svg` file.
I've included those files in this repository: `x-tree.svg` and `x-arc.svg`.
These are converted to .eps files by the `Makefile` in the project parent
directory.

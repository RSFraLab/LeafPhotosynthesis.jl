Vendored browser libraries of the LeafPhotosynthesis.jl manual
================================================================

These files are used only by the interactive charts of the documentation
(docs/src/public/leaf/). They add no dependency to the LeafPhotosynthesis Julia
package and are served from the site itself, so the charts do not need a CDN.

Chart.js 4.4.8
--------------
chart.umd.js: Chart.js 4.4.8 UMD bundle, vendored unchanged from
https://cdn.jsdelivr.net/npm/chart.js@4.4.8/dist/chart.umd.js
Upstream: https://github.com/chartjs/Chart.js/tree/v4.4.8
License: MIT, see Chart.js-LICENSE.txt. The bundled @kurkle/color 0.3.2
dependency is also MIT licensed, see kurkle-color-LICENSE.txt.

KaTeX 0.16.47
-------------
katex/: KaTeX 0.16.47 CSS, core renderer, auto-render extension and WOFF2 math
fonts, vendored unchanged from the katex@0.16.47 npm distribution. Used for the
equations inside the Flexas et al. (2002) explorer.
Upstream: https://github.com/KaTeX/KaTeX
License: MIT, see katex/KaTeX-LICENSE.txt.

The licence files were renamed from .md to .txt (content unchanged) so that the
documentation build does not treat them as pages.

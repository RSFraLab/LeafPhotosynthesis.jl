# Offline documentation dependencies

## Chart.js 4.4.8

Chart.js 4.4.8 UMD bundle (`chart.umd.js`), vendored unchanged from:
https://cdn.jsdelivr.net/npm/chart.js@4.4.8/dist/chart.umd.js

Upstream: https://github.com/chartjs/Chart.js/tree/v4.4.8
License: MIT; see `Chart.js-LICENSE.md`. The bundled `@kurkle/color` dependency is
also MIT licensed; see `kurkle-color-LICENSE.md` (version 0.3.2).

Chart integration and accessibility guidance:
https://www.chartjs.org/docs/latest/getting-started/integration.html
https://www.chartjs.org/docs/latest/general/accessibility.html

## KaTeX 0.16.47

KaTeX 0.16.47 CSS, core renderer, auto-render extension, `math/tex` script extension,
and WOFF2 math fonts (`katex/`), vendored unchanged from the `katex@0.16.47` npm
distribution for offline equation rendering on GitHub Pages (`docs/_includes/head-custom.html`).

Upstream: https://github.com/KaTeX/KaTeX
License: MIT; see `katex/KaTeX-LICENSE.md`.

These pinned browser assets are used strictly by the documentation and add zero
dependencies to the `LeafPhotosynthesis` Julia package.

# Offline plotting dependency

Chart.js 4.4.8 UMD bundle, vendored unchanged from:
https://cdn.jsdelivr.net/npm/chart.js@4.4.8/dist/chart.umd.js

Upstream: https://github.com/chartjs/Chart.js/tree/v4.4.8
License: MIT; see Chart.js-LICENSE.md. The bundled @kurkle/color dependency is
also MIT licensed; see kurkle-color-LICENSE.md (version 0.3.2).

This pinned browser dependency is only used by the documentation. It does not add
any dependency to the Julia model. The source-map reference in the upstream bundle
is retained; a source map is not required to view the plots.

Chart integration and accessibility guidance:
https://www.chartjs.org/docs/latest/getting-started/integration.html
https://www.chartjs.org/docs/latest/general/accessibility.html

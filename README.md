# Vietoris-Rips complex builder

Small library for building Vietoris-Rips Complex using incremental algorithm as described in the following paper:
[*"Fast Construction of the Vietoris-Rips Complex"*](http://www.cs.sandia.gov/CSRI/Workshops/2009/CAT/presentations/zomorodian.pdf)
by Afra Zomorodian.

Construction of Vietoris-Rips Complex is much faster than Cech Complex so it has more practical applications.

## Installation

```bash
npm install vr-complex
```

or in browser

```bash
bower install vr-complex
```

## Usage

Building complex (Betti numbers)
```js
var VR = require('vr-complex');
var vr = new VR();

var cells = [
  [0, 0], // [x, y]
  [1, 1],
  // ...
];
var maxK = 3; // maximum size of simplex (0 - point, 1 - edge, 2 - triangle, ...)
var R = 10; // radius

var simplices = vr.complex(cells, maxK, R);
/*
[
  [[0], [1], ... ], // 0-simplex (points)
  [[0,1], [0,5], ...], // 1-simplex (edges)
  [[0,1,5], [0,1,9], ...], // 2-simplex (triangles)
  ... // 3-simplex (tetrahedrons)
]
*/

```

Computing homology
```js
var simplices = vr.complex(cells, maxK, R); // compute simplices
var homology = vr.homology(3); // compute first 3 Betti numbers
```

## Test

To test, install `mocha` globally and run following command:
```bash
npm test
```

## Future releases
- benchmark

/**
 * Vietrois-Rips complex builder
 */
function VR() {
  this.S;
}

/* top level functions */

/**
 * Build VR complex
 *
 * @param {Array.<Array>} data
 * @param {number} maxK
 * @param {number} R
 * @returns {Array.<Array>}
 */
VR.prototype.complex = function(data, maxK, R) {
  var neighbors = this.neighborhood(data, R);
  var vr = this.incrementalVR(neighbors, maxK);
  var l = vr.length;

  // Group simplices by dimension:
  // S[0] - points
  // S[1] - edges
  // S[2] - triangles
  // ...
  this.S = [];
  for (var i = 0; i <= maxK; i++) {
    this.S.push([]);
  }

  while (l--) {
    var simplex = vr[l];
    this.S[simplex.length - 1].push(simplex);
  }

  return this.S;
};

/**
 * Compute homology (betti numbers)
 *
 * @param {number} maxK
 * @returns {Array.<number>}
 */
VR.prototype.homology = function(maxK) {
  if (typeof maxK === 'undefined') {
    maxK = this.S.length - 1;
  }

  var bd = [];
  var betti = [];
  var bounds;

  for (var i = 0; i <= maxK; i++) {
    if (!this.S[i].length) {
      break;
    }

    if (i === 0) {
      bounds = [this.zeros(this.S[i].length)];
    } else {
      bounds = this.boundaryMatrix(this.S[i - 1], this.S[i]);
      betti.push(this.bettiNumber(bd[i - 1], bounds));
    }
    bd.push(bounds);
  }

  return betti;
};

/* core of VR algorithm */

/**
 * Incrementsl Vietrois-Rips algorithm
 *
 * @param {Array.<Array>} G - neighborhood graph
 * @param {number} maxK
 * @returns {Array.<Array>}
 */
VR.prototype.incrementalVR = function(G, k) {
  var V = [];
  for (var u = 0, l = G.length; u < l; u++) {
    this.addCofaces(G, k, [u], G[u], V);
  }
  return V;
};

/**
 * Recurrent part of incremental algorithm
 *
 * @param {Array} G - neighborhood graph
 * @param {number} k - maximum simplex dimension
 * @param {Array} t - simplex elements
 * @param {Array} N - array of vertex neighbors
 * @param {Array} V - output, list of simplices
 * @returns {undefined}
 */
VR.prototype.addCofaces = function(G, k, t, N, V) {
  V.push(t);
  if (t.length > k) {
    return;
  } else {
    for (var i = 0; i < N.length; i++) {
      var v = N[i];
      var s = [v].concat(t);
      var M = this.intersection(N, G[v]);
      this.addCofaces(G, k, s, M, V);
    }
  }
};

/**
 * Building neighborhood graph for given set of points
 *
 * @param {Array.<Array>} S - set of points
 * @param (number> e - distance threshold
 * @returns {Array.<Array>}
 */
VR.prototype.neighborhood = function(S, e) {
  var graph = [];
  var i = S.length;

  while (i-- > 0) {
    var neighbors = [];
    var j = i;
    while (j--) {
      if (this.distance(S[i], S[j]) <= e) {
        neighbors.unshift(j);
      }
    }
    graph.unshift(neighbors);
  }

  return graph;
};

/**
 * Euclidean distance
 * Only for low-dimensional spaces ( < 20)
 *
 * @param {number} maxK
 * @returns {Array.<number>}
 */
VR.prototype.distance = function(p, q) {
  var sum = 0;
  var i = p.length;
  while (i--) {
    sum += (p[i] - q[i]) * (p[i] - q[i]);
  }
  return Math.sqrt(sum);
};

/**
 * Intersection of two sets
 * Important: both sets must be sorted!
 *
 * @param {Array.<number>} a
 * @param {Array.<number>} b
 * @returns {Array.<number>}
 */
VR.prototype.intersection = function(a, b) {
  var ai = 0;
  var bi = 0;
  var result = [];

  while (ai < a.length && bi < b.length) {
    if (a[ai] < b[bi]) {
      ai++;
    } else if (a[ai] > b[bi]) {
      bi++;
    } else {
      result.push(a[ai]);
      ai++;
      bi++;
    }
  }

  return result;
};

/**
 * Compute k-th Betti number
 * k-th betti number is defined as a difference of
 * image and kernel dimension of two boundary matrices
 *
 * @param {Array.<Array>} k-th boundary matrix
 * @param {Array.<Array>} (k+1)-th boundary matrix
 * @returns {Array.<number>}
 */
VR.prototype.bettiNumber = function(dK, dKplus1) {
  var A = dK.slice();
  var B = dKplus1.slice();
  this.simultaneousReduce(A, B);
  this.finishRowReducing(B);

  var dimKChains = A[0].length;
  var kernelDim = dimKChains - this.numPivotCols(A);
  var imageDim = this.numPivotRows(B);

  return kernelDim - imageDim;
};

/**
 * Build boundary matrix
 *
 * Boundary matrix is a two dimensional array,
 * where columns labels are (k+1)-simplices and row labels are k-simplices.
 * Each value represents sign of hat() operator on a given (k+1) simplex
 * @example
 * hat([1,2,3],1) = -1 * [1,3]
 * -> Boundary[[1,3]],[[1,2,3]] will have value -1
 *
 * @param {Array.<Array>} A - input array
 * @param {number} i
 * @param {number} j
 * @returns {Array.<Array>}
 */
VR.prototype.boundaryMatrix = function(rows, cols) {
  var numRows = rows.length;
  var numCols = cols.length;
  var i;

  // build hash table for faster computation
  var hash = {};
  for (i = 0 ; i < numRows; i++) {
    hash[rows[i]] = i;
  }

  var A = this.zeros(numRows, numCols);

  for (var colId = 0; colId < numCols; colId++) {
    var v = cols[colId];
    for (i = 0, l = v.length; i < l; i++) {
      var hatted = this.hatOperator(v, i);
      var rowId = -1;
      if (hash.hasOwnProperty(hatted[1])) {
        rowId = hash[hatted[1]];
      }

      A[rowId][colId] = (rowId >= 0) ? hatted[0] : 0;
    }
  }

  return A;
};

/* linear algebra utils */

/**
 * Hat operator
 * Removes element of given id from simplex.
 * Main purpose is to divide simplex into simplices of lower dimension
 *
 * @example
 * triangle [1,2,3] can be split into 3 2-simplices:
 * [2,3], -[1,3], [1,2]
 * (-1) in front of simplex reverse it: -[1,3] == [3,1]
 *
 * @param {Array.<number>} simplex
 * @param {number} i - id of element to remove
 * @returns {Array.<number>}
 */
VR.prototype.hatOperator = function(simplex, i) {
  var s = simplex.slice();
  s.splice(i, 1);
  return [Math.pow(-1, i), s];
};

/**
 * Reduce two boundary matrices into echelon form simultaneously
 *
 * @param {Array.<Array>} A
 * @param {Array.<Array>} B
 * @returns {Array.<Array>}
 */
VR.prototype.simultaneousReduce = function(A, B) {
  if (A[0].length !== B.length) {
    throw new Error('Matrices have the wrong shape.');
  }

  var numRows = A.length;
  var numCols = A[0].length;
  var i = 0;
  var j = 0;

  while (true) {
    if (i >= numRows || j >= numCols) {
      break;
    }

    if (A[i][j] === 0) {
      var nonzeroCol = j;

      while (nonzeroCol < numCols && A[i][nonzeroCol] === 0) {
        nonzeroCol += 1;
      }

      if (nonzeroCol == numCols) {
        i++;
        continue;
      }

      this.colSwap(A, j, nonzeroCol);
      this.rowSwap(B, j, nonzeroCol);
    }

    var pivot = A[i][j];

    this.scaleCol(A, j, 1.0 / pivot);
    this.scaleRow(B, j, 1.0 / pivot);

    for (var otherCol = 0; otherCol < numCols; otherCol++) {
      if (otherCol === j) {
        continue;
      }

      if (A[i][otherCol] !== 0) {
        scaleAmt = -A[i][otherCol];
        this.colCombine(A, otherCol, j, scaleAmt);
        this.rowCombine(B, j, otherCol, -scaleAmt);
      }
    }

    i++;
    j++;
  }

  return [A, B];
};

/**
 * After simultaneous reduction, small part of one array still has to be cleaned up
 *
 * @param {Array.<Array>} B - boundary matrix
 * @returns {Array.<Array>}
 */
VR.prototype.finishRowReducing = function(B) {
  var numRows = B.length;
  var numCols = B[0].length;

  var i = 0;
  var j = 0;
  var nonzeroRow;

  while (true) {
    if (i >= numRows || j >= numCols) {
      break;
    }

    if (B[i][j] === 0) {
      nonzeroRow = i;
      while (nonzeroRow < numRows && B[nonzeroRow][j] === 0) {
        nonzeroRow += 1;
      }

      if (nonzeroRow === numRows) {
        j++;
        continue;
      }
      this.rowSwap(B, i, nonzeroRow);
    }
    pivot = B[i][j];
    this.scaleRow(B, i, 1.0 / pivot);

    for (var otherRow = 0; otherRow < numRows; otherRow++) {
      if (otherRow == i) {
        continue;
      }

      if (B[otherRow][j] !== 0) {
        scaleAmt = -B[otherRow][j];
        this.rowCombine(B, otherRow, i, scaleAmt);
      }
    }

    i++;
    j++;
  }

  return B;
};

/**
 * Swap two rows
 *
 * @param {Array.<Array>} A - input array
 * @param {number} i
 * @param {number} j
 * @returns {Array.<Array>}
 */
VR.prototype.rowSwap = function(A, i, j) {
  var temp = A[i];
  A[i] = A[j];
  A[j] = temp;
  return A;
};

/**
 * Swap two columns
 *
 * @param {Array.<Array>} A - input array
 * @param {number} i
 * @param {number} j
 * @returns {Array.<Array>}
 */
VR.prototype.colSwap = function(A, i, j) {
  var temp;
  for (var k = 0, len = A.length; k < len; k++) {
    temp = A[k][i];
    A[k][i] = A[k][j];
    A[k][j] = temp;
  }
  return A;
};

/**
 * Multiply row by a value
 *
 * @param {Array.<Array>} A - input array
 * @param {number} colId
 * @param {number} value
 * @returns {Array.<Array>}
 */
VR.prototype.scaleRow = function(A, rowId, multiplier) {
  for (var j = 0, len = A[rowId].length; j < len; j++) {
    A[rowId][j] *= multiplier;
  }
  return A;
};

/**
 * Multiply column by a value
 *
 * @param {Array.<Array>} A - input array
 * @param {number} colId
 * @param {number} value
 * @returns {Array.<Array>}
 */
VR.prototype.scaleCol = function(A, colId, multiplier) {
  for (var j = 0, len = A.length; j < len; j++) {
    A[j][colId] *= multiplier;
  }
  return A;
};

/**
 * Multiply row by a value and add to another row
 *
 * @param {Array.<Array>} A - input array
 * @param {number} addTo - id of the destination row
 * @param {number} scaleRow - id of source row
 * @param {number} multiplier
 * @returns {Array.<Array>}
 */
VR.prototype.rowCombine = function(A, addTo, scaleRow, multiplier) {
  for (var j = 0, len = A[addTo].length; j < len; j++) {
    A[addTo][j] += multiplier * A[scaleRow][j];
  }
  return A;
};

/**
 * Multiply column by a value and add to another column
 *
 * @param {Array.<Array>} A - input array
 * @param {number} addTo - id of the destination column
 * @param {number} scaleCol - id of source column
 * @param {number} multiplier
 * @returns {Array.<Array>}
 */
VR.prototype.colCombine = function(A, addTo, scaleCol, multiplier) {
  for (var j = 0, len = A.length; j < len; j++) {
    A[j][addTo] += multiplier * A[j][scaleCol];
  }
  return A;
};

/**
 * Count non-zero rows
 *
 * @param {Array.<Array>} A
 * @returns {number}
 */
VR.prototype.numPivotRows = function(A) {
  var z = 0;
  for (var i = 0, len1 = A.length; i < len1; i++) {
    for (var j = 0, len2 = A[i].length; j < len2; j++) {
      if (A[i][j] !== 0) {
        z++;
        break;
      }
    }
  }
  return z;
};

/**
 * Count non-zero columns
 *
 * @param {Array.<Array>} A
 * @returns {number}
 */
VR.prototype.numPivotCols = function(A) {
  var z = 0;
  for (var j = 0, len1 = A[0].length; j < len1; j++) {
    for (var i = 0, len2 = A.length; i < len2; i++) {
      if (A[i][j] !== 0) {
        z++;
        break;
      }
    }
  }
  return z;
};

/* utils */

/**
 * Create array filled with zeros
 *
 * @param {number} numRows
 * @param {nummber} numCols
 * @returns {Array}
 */
VR.prototype.zeros = function(numRows, numCols) {
  var A = [];

  for (var i = 0; i < numRows; i++) {
    if (typeof numCols !== 'undefined') {
      A.push(this.zeros(numCols));
    } else {
      A.push(0);
    }
  }
  return A;
};

// publish module
if (typeof module !== 'undefined' && typeof module.exports !== 'undefined') {
  module.exports = VR;
}

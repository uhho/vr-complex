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

  this.S = [];
  for (var i = 0; i <= maxK; i++) {
    this.S.push([]);
  }

  var l = vr.length;

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
 * Incremental VR-Complex algorithm
 * Based on "Fast Construction of the Vietoris-Rips Complex"
 * by Afra Zomorodian
 *
 * @param  {[type]} G [description]
 * @param  {[type]} k [description]
 * @return {[type]}   [description]
 */
VR.prototype.incrementalVR = function(G, k) {
  var V = [];
  for (var u = 0, l = G.length; u < l; u++) {
    this.addCofaces(G, k, [u], G[u], V);
  }
  return V;
};

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

VR.prototype.neighborhood = function(S, e) {
  var graph = [];
  var i = S.length;
  var d;

  while (i-- > 0) {
    var neighbors = [];
    var j = i;
    while (j--) {
      d = this.distance(S[i], S[j]);
      if (d <= e) {
        neighbors.unshift(j);
      }
    }
    graph.unshift(neighbors);
  }

  return graph;
};

VR.prototype.distance = function(p, q) {
  var sum = 0;
  var i = p.length;
  while (i--) {
    sum += (p[i] - q[i]) * (p[i] - q[i]);
  }
  return Math.sqrt(sum);
};

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

VR.prototype.hatOperator = function(simplex, i) {
  var s = simplex.slice();
  s.splice(i, 1);
  return [Math.pow(-1, i), s];
};

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

VR.prototype.rowSwap = function(A, i, j) {
  var temp = A[i];
  A[i] = A[j];
  A[j] = temp;
  return A;
};

VR.prototype.colSwap = function(A, i, j) {
  var temp;
  for (var k = 0, len = A.length; k < len; k++) {
    temp = A[k][i];
    A[k][i] = A[k][j];
    A[k][j] = temp;
  }
  return A;
};

VR.prototype.scaleRow = function(A, i, c) {
  for (var j = 0, len = A[i].length; j < len; j++) {
    A[i][j] *= c;
  }
  return A;
};

VR.prototype.scaleCol = function(A, i, c) {
  for (var j = 0, len = A.length; j < len; j++) {
    A[j][i] *= c;
  }
  return A;
};

VR.prototype.rowCombine = function(A, addTo, scaleRow, scaleAmt) {
  for (var j = 0, len = A[addTo].length; j < len; j++) {
    A[addTo][j] += scaleAmt * A[scaleRow][j];
  }
  return A;
};

VR.prototype.colCombine = function(A, addTo, scaleCol, scaleAmt) {
  for (var j = 0, len = A.length; j < len; j++) {
    A[j][addTo] += scaleAmt * A[j][scaleCol];
  }
  return A;
};

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

VR.prototype.getArrayId = function(arr, el) {
  var len = arr.length;
  for (var i = 0; i < len; i++) {
    if (arr[i].join() === el.join()) {
      return i;
    }
  }
  return -1;
};

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

if (typeof module !== 'undefined' && typeof module.exports !== 'undefined') {
  module.exports = VR;
}

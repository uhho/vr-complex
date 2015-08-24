var VR = require('../lib/VR.js');

describe('VR', function() {
  describe('complex', function() {
    it('should create correct VR complex', function(done) {
      var data = [
        [0, 0],
        [0, 1],
        [1, 0],
        [1, 1]
      ];
      var vr = new VR();
      var S = vr.complex(data, 3, 3);
      S[3].should.eql([[0, 1, 2, 3]]);
      S.length.should.eql(4);
      done();
    });
  });

  describe('bettiNumbers', function() {
    it('should compute correct betti numbers - no hole', function(done) {
      var data = [
        [-1, 0],
        [1, 0],
        [0, 2],
        [0, -2]
      ];
      var vr = new VR();
      var S = vr.complex(data, 2, 3);
      var betti = vr.homology(2);
      betti.should.eql([1, 0]);
      done();
    });

    it('should compute correct betti numbers - hole', function(done) {
      var data = [
        [6.06, 2.01], [5.10, 0.10], [4.09, 2.02],
        [3.10, 4.02],
        [2.01, 6.06], [2.01, 8.08], [4.05, 7.07],
        [6.06, 7.07],
        [8, 6], [9, 6], [11, 7], [9, 4], [7, 4]
      ];
      var vr = new VR();
      var S = vr.complex(data, 2, 3);
      var betti = vr.homology(2);
      betti.should.eql([1, 1]);
      done();
    });

  });

  describe('distance', function() {
    it('should compute distance between two points', function(done) {
      var vr = new VR();
      vr.distance([0 ,0], [1, 0]).should.eql(1);
      vr.distance([0 ,-1], [0, 1]).should.eql(2);
      vr.distance([-1 ,-1], [1, 1]).should.eql(2 * Math.sqrt(2));
      done();
    });
  });

  describe('intersection', function() {
    it('should return intersection of two sets', function(done) {
      var vr = new VR();
      vr.intersection([], []).should.eql([]);
      vr.intersection([1, 2, 3], [1, 2, 3]).should.eql([1, 2, 3]);
      vr.intersection([0, 1, 2, 3, 4], [0, 2, 4]).should.eql([0, 2, 4]);
      vr.intersection([0, 1, 1], [0, 1, 2, 2]).should.eql([0, 1]);
      done();
    });
  });

  describe('bettiNumber', function() {
    it('should compute k-th betti number');
  });

  describe('boundaryMatrix', function() {
    it('should comput boundary matrix', function(done) {
      var rows = [[1, 2], [1, 4], [2, 3], [2, 4], [3, 4]];
      var cols = [[1, 2, 4], [2, 3, 4]];
      var vr = new VR();
      var bd = vr.boundaryMatrix(rows, cols);
      bd.should.eql([
        [1, 0],
        [-1, 0],
        [0, 1],
        [1, -1],
        [0, 1]
      ]);
      done();
    });
  });

  describe('hatOperator', function() {
    it('should remove element from vector', function(done) {
      var vr = new VR();
      vr.hatOperator([1, 2, 3, 4, 5], 2).should.eql([1, [1, 2, 4, 5]]);
      done();
    });
  });

  describe('simultaneousReduce', function() {
    it('should reduce two boundary matrices into echelon form');
  });

  describe('finishRowReducing', function() {
    it('should finish matrix reduction');
  });

  describe('rowSwap', function() {
    it('should swap two rows', function(done) {
      var vr = new VR();
      var a = [
        [2, 1, 3],
        [0, 0, 0],
        [1, 3, 2]
      ];
      vr.rowSwap(a, 0, 2).should.eql([
        [1, 3, 2],
        [0, 0, 0],
        [2, 1, 3]
      ]);
      done();
    });
  });

  describe('colSwap', function() {
    it('should swap two columns', function(done) {
      var vr = new VR();
      var a = [
        [1, 0, 2],
        [3, 0, 1],
        [2, 0, 3]
      ];
      vr.colSwap(a, 0, 2).should.eql([
        [2, 0, 1],
        [1, 0, 3],
        [3, 0, 2]
      ]);
      done();
    });
  });

  describe('scaleRow', function() {
    it('should swap two rows', function(done) {
      var vr = new VR();
      var a = [
        [2, 1, 3],
        [0, 0, 0],
        [1, 3, 2]
      ];
      vr.scaleRow(a, 2, 2).should.eql([
        [2, 1, 3],
        [0, 0, 0],
        [2, 6, 4]
      ]);
      done();
    });
  });

  describe('scaleCol', function() {
    it('should swap two columns', function(done) {
      var vr = new VR();
      var a = [
        [1, 0, 2],
        [3, 0, 1],
        [2, 0, 3]
      ];
      vr.scaleCol(a, 2, 2).should.eql([
        [1, 0, 4],
        [3, 0, 2],
        [2, 0, 6]
      ]);
      done();
    });
  });

  describe('rowCombine', function() {
    it('should find number of non-zero rows', function(done) {
      var vr = new VR();
      var a = [
        [2, 1, 3],
        [0, 0, 0],
        [1, 3, 2]
      ];
      var result = vr.rowCombine(a, 0, 2, 2).should.eql([
        [4, 7, 7],
        [0, 0, 0],
        [1, 3, 2]
      ]);
      done();
    });
  });

  describe('colCombine', function() {
    it('should find number of non-zero rows', function(done) {
      var vr = new VR();
      var a = [
        [1, 0, 2],
        [3, 0, 1],
        [2, 0, 3]
      ];
      vr.colCombine(a, 0, 2, 2).should.eql([
        [5, 0, 2],
        [5, 0, 1],
        [8, 0, 3]
      ]);
      done();
    });
  });

  describe('numPivotRows', function() {
    it('should find number of non-zero rows', function(done) {
      var vr = new VR();
      var a = [
        [2, 1, 3],
        [0, 0, 0],
        [1, 3, 2]
      ];
      vr.numPivotRows(a).should.eql(2);
      done();
    });
  });

  describe('numPivotCols', function() {
    it('should find number of non-zero cols', function(done) {
      var vr = new VR();
      var a = [
        [1, 0, 2],
        [3, 0, 1],
        [2, 0, 3]
      ];
      vr.numPivotCols(a).should.eql(2);
      done();
    });
  });

  describe('zeros', function() {
    it('should create 0-array', function(done) {
      var vr = new VR();
      vr.zeros(2).should.eql([0, 0]);
      vr.zeros(2, 3).should.eql([
        [0, 0, 0],
        [0, 0, 0]
      ]);
      done();
    });
  });

});

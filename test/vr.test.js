var VR = require('../VR.js');

describe('Cech', function() {
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
    it('should compute correct betti numbers', function(done) {
      var data = [
        [-1, 0],
        [1, 0],
        [0, 2],
        [0, -2]
      ];
      var vr = new VR();
      var S = vr.complex(data, 2, 3);
      var betti = vr.homology(S);
      betti.should.eql([1, 0]);
      done();
    });
    
    it('should compute correct betti numbers', function(done) {
      var data = [
        [6.06, 2.01], [5.10, 0.10], [4.09, 2.02],
        [3.10, 4.02],
        [2.01, 6.06], [2.01, 8.08], [4.05, 7.07],
        [6.06, 7.07],
        [8, 6], [9, 6], [11, 7], [9, 4], [7, 4]
      ];
      var vr = new VR();
      var S = vr.complex(data, 2, 3);
      var betti = vr.homology(S);
      betti.should.eql([1, 1]);
      done();
    });
    
  });
  
  describe('boundaryMatrix', function() {
    it('should create boundary matrix');
  });
  
   describe('hatOperator', function() {
    it('should remove element from vector', function(done) {
      vr.hatOperator([1, 2, 3, 4, 5], 2).should.eql([1, 2, 4, 5]);
      done();
    });
  });

  
});

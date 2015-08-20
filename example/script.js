var $radiusLabel = $('#radiusLabel');
var $radius = $('#radius');
var $betti0 = $('#betti0');
var $betti1 = $('#betti1');
var $shadow = $('#shadow');
var $fps = $('#fps');
var $vertices = $('#vertices');
var $edges = $('#edges');
var $triangles = $('#triangles');
var $tetrahedrons = $('#tetrahedrons');

$radius.on('input', function() {
  var R = $(this).val();
  $radiusLabel.text(R);
});

var fps = 0;
var now;
var lastUpdate = (new Date()) * 1;
var fpsFilter = 50;

var canvas = document.getElementById('canvas');
var ctx = canvas.getContext('2d');
var width = canvas.width;
var height = canvas.height;

/******************************************************************************/

// generate data
var vr = new VR();
var cellR = 90;
var a = generateData(200, 200, 160, 11, cellR);
var b = generateData(200, 200, 120, 9, cellR);
var c = generateData(440, 200, 140, 9, cellR);
var d = generateData(440, 200, 100, 7, cellR);
var cells = a.concat(b, c, d);

// update complex
setInterval(function() {
  var R = parseInt($radius.val());
  var l = cells.length;
  var cell;
  var variance = 5;
  while (l--) {
    cell = cells[l];

    cell[0] += (Math.random() * variance) - variance / 2;
    cell[1] += (Math.random() * variance) - variance / 2;

    cell[0] = Math.min(width, Math.max(0, cell[0]));
    cell[1] = Math.min(height, Math.max(0, cell[1]));
  }

  draw(cells, R);
}, 50);

// update homology
setInterval(function() {
  var betti = vr.homology(2);
  $betti0.text(betti[0]);
  $betti1.text(betti[1]);
}, 100);

// updated FPS
setInterval(function() {
  $fps.text(fps.toFixed(1) + 'fps');
}, 1000);

/******************************************************************************/

function draw(cells, R) {
  var S = vr.complex(cells, 3, R);

  $vertices.text(S[0].length);
  $edges.text(S[1].length);
  $triangles.text(S[2].length);
  $tetrahedrons.text(S[3].length);

  ctx.clearRect(0, 0, width, height);
  // draw the circle
  drawCircles(cells, R);
  drawPaths(S[1], true, false);
  drawPaths(S[2], false, true);

  var thisFrameFPS = 1000 / ((now = new Date()) - lastUpdate);
  if (now != lastUpdate) {
    fps += (thisFrameFPS - fps) / fpsFilter;
    lastUpdate = now;
  }
}

function drawCircles(cells, R) {
  var l = cells.length;
  for (var i = 0; i < l; i++) {
    ctx.globalAlpha = 0.02;
    ctx.fillStyle = '#000000';
    ctx.strokeStyle = '#000000';
    ctx.beginPath();
    ctx.arc(cells[i][0], cells[i][1], R / 2, 0, 2 * Math.PI, false);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    ctx.globalAlpha = 0.5;
    ctx.fillStyle = '#ED7B3D';
    ctx.strokeStyle = '#ED7B3D';
    ctx.beginPath();
    ctx.arc(cells[i][0], cells[i][1], 3, 0, 2 * Math.PI, false);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
  }
}

function drawPaths(S, stroke, fill) {
  if ($shadow.hasClass('active')) {
    ctx.globalAlpha = 1;
    ctx.fillStyle = '#555555';
    ctx.strokeStyle = '#555555';
  } else {
    ctx.globalAlpha = 0.1;
    ctx.fillStyle = '#74197D';
    ctx.strokeStyle = '#74197D';
  }

  ctx.lineWidth = 2;
  l = S.length;
  for (var i = 0; i < l; i++) {
    var simplex = S[i];
    ctx.beginPath();
    ctx.moveTo(cells[simplex[0]][0], cells[simplex[0]][1]);
    for (var j = 1; j < simplex.length; j++) {
      ctx.lineTo(cells[simplex[j]][0], cells[simplex[j]][1]);
    }
    ctx.closePath();
    if (stroke) {
      ctx.stroke();
    }
    if (fill) {
      ctx.fill();
    }
  }
}

function generateData(h, k, R, len, r) {
  var cells = [];
  var step = (2 * Math.PI / len) + Math.random() * 0.1;
  var twoPI = 2 * Math.PI;

  for (var theta = 0; theta < twoPI; theta += step) {
    var x = h + R * Math.cos(theta);
    x += (Math.random() * 20) - 10;
    var y = k - R * Math.sin(theta);
    y += (Math.random() * 20) - 10;

    cells.push([x, y]);
  }

  return cells;
}

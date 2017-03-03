var variantList = [];
var geneNames = [];
var maxCov = 0;
var maxQual = 0;
var geneMode = 0;
var windowSize = [0, 0];
var zoomSize = [0, 0];
var dragOn = false;
var covThreshold = [0, 0];
var displayCoverage = false;
var freqThreshold = [0, 0];
var displayFrequency = false;
var qualThreshold = [0, 0];
var displayQuality = false;
var biasThreshold = [0, 0];
var displayBias = false;
var geneEndPoints = [];
var geneChrs = [];

function readVCF(evt) {
	var file = evt.target.files[0];
	var reader = new FileReader();
	reader.onload = function() {
		variantList = [];
		var lines = reader.result.split(/[\r\n]+/g);
		for (i = 1; i < lines.length; i ++) {
			var pieces = lines[i].split('\t');
			if (geneNames.indexOf(pieces[6]) < 0) {
				geneNames.push(pieces[6]);
			}
			var c = "";
			if (pieces[0].length == 5) {
				c = pieces[0].charAt(3) + pieces[0].charAt(4);
			}
			else {
				c = pieces[0].charAt(3);
			}
			variantList.push( { chr:parseInt(c), pos:parseInt(pieces[1]), refN:pieces[2], varN:pieces[3], freq:parseFloat(pieces[4]), qual:parseFloat(pieces[5]), geneID:pieces[6], cov:parseInt(pieces[7]), strandBias:parseFloat(pieces[9]), name:pieces[10], phred:parseFloat(pieces[11]), cadd:parseInt(pieces[12]) } );

			if (parseInt(pieces[7]) > maxCov) {
				maxCov = parseInt(pieces[7]);
			}		
			if (parseInt(pieces[5]) > maxQual) {
				maxQual = parseInt(pieces[5]);
			} 
		}

		document.getElementById('covRangeLower').max = maxCov;
		document.getElementById('covRangeUpper').max = maxCov;
		document.getElementById('covRangeUpper').value = maxCov;
		document.getElementById('covTextUpper').value = maxCov;
		covThreshold[1] = maxCov;
		document.getElementById('freqRangeUpper').value = 100;
		document.getElementById('freqTextUpper').value = 100;
		freqThreshold[1] = 100;
		document.getElementById('qualRangeLower').max = maxQual;
		document.getElementById('qualRangeUpper').max = maxQual;
		document.getElementById('qualRangeUpper').value = maxQual;
		document.getElementById('qualTextUpper').value = maxQual;
		qualThreshold[1] = maxQual;
		document.getElementById('biasRangeUpper').value = 1;
		document.getElementById('biasTextUpper').value = 1;
		biasThreshold[1] = 1;

		var nuc = document.getElementById("nucPanel");
		var nx = nuc.getContext("2d");
		nx.strokeStyle = "#000000";
		nx.beginPath();
		nx.moveTo(0, 90);
		nx.lineTo(20, 90);
		nx.moveTo(0, 105);
		nx.lineTo(20, 105);
		nx.moveTo(0, 120);
		nx.lineTo(20, 120);
		nx.moveTo(0, 135);
		nx.lineTo(20, 135);
		nx.stroke();
		nx.font = "10px Georgia";
		nx.fillStyle = "#FF0000";
		nx.fillText("G", 5, 100);
		nx.fillStyle = "#AFFF00";
		nx.fillText("C", 5, 115);
		nx.fillStyle = "#0000FF";
		nx.fillText("T", 5, 130);
		nx.fillStyle = "#00FFFF";
		nx.fillText("A", 5, 145);


		geneNames.sort();
		var s = document.getElementById('geneName');
		for (i = 0; i < geneNames.length; i ++) {
			var option = document.createElement("option");
			option.text = geneNames[i];
			s.add(option);
			geneEndPoints.push([10000000000, 0]);
		}
		s.options[0].selected;

		for (i = 0; i < variantList.length; i++) {
			var geneNum = geneNames.indexOf(variantList[i].geneID);
			if (variantList[i].pos < geneEndPoints[geneNum][0]) {
				geneEndPoints[geneNum][0] = variantList[i].pos;
			}
			if (variantList[i].pos > geneEndPoints[geneNum][1]) {
				geneEndPoints[geneNum][1] = variantList[i].pos;
			}
			geneChrs[geneNum] = variantList[i].chr;
		}
		for (i = 0; i < geneEndPoints.length; i ++) {
			geneEndPoints[i][0] -= 1000;
			geneEndPoints[i][1] += 1000;
		}

		document.getElementById('ucsc').style.visibility = "hidden";
		document.getElementById('ucsc').src="https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position=chr" + geneChrs[0] + ":" + geneEndPoints[0][0] + "-" + geneEndPoints[0][1];
		document.getElementById('chrNum').value = "chr" + geneChrs[0];
		document.getElementById('start').value = geneEndPoints[0][0];
		document.getElementById('end').value = geneEndPoints[0][1];
		windowSize[0] = geneEndPoints[0][0];
		windowSize[1] = geneEndPoints[0][1];
		zoomSize[0] = geneEndPoints[0][0];
		zoomSize[1] = geneEndPoints[0][1];
	};
	reader.readAsText(file);
}

function draw() {

	var canvas = document.getElementById("mainTrack");
	var ctx = canvas.getContext("2d");
	ctx.clearRect(0, 0, 660, 150);
	ctx.font = "10px Georgia";

	var nucView = document.getElementById("nucPanel");
	var nx = nucView.getContext("2d");
	nx.clearRect(0, 0, 20, 90);
	nx.strokeStyle = "#000000";

	var cad = document.getElementById("cadd");
	var cdx = cad.getContext("2d");
	cdx.clearRect(0, 0, 660, 25);
	cdx.fillStyle = "#FFFFFF";
	cdx.strokeStyle = "#000000";

	var prot = document.getElementById("polyphen");
	var px = prot.getContext("2d");
	px.clearRect(0, 0, 660, 25);
	px.fillStyle = "#FFFFFF";
	px.strokeStyle = "#000000";

	ctx.strokeStyle = "#000000";
	ctx.beginPath();
	ctx.moveTo(0, 90);
	ctx.lineTo(660, 90);
	ctx.moveTo(0, 105);
	ctx.lineTo(660, 105);
	ctx.moveTo(0, 120);
	ctx.lineTo(660, 120);
	ctx.moveTo(0, 135);
	ctx.lineTo(660, 135);
	ctx.stroke();

	var x = 0;
	var y = 0;
	for (i = 0; i < variantList.length; i ++) {
		if (variantList[i].geneID == geneNames[geneMode] && variantList[i].pos >= zoomSize[0] && variantList[i].pos <= zoomSize[1]) {
			if (variantList[i].cov >= covThreshold[0] && variantList[i].cov <= covThreshold[1] && variantList[i].freq >= freqThreshold[0] && variantList[i].freq <= freqThreshold[1] && variantList[i].qual >= qualThreshold[0] && variantList[i].qual <= qualThreshold[1] && variantList[i].strandBias >= biasThreshold[0] && variantList[i].strandBias <= biasThreshold[1]) {
				var posX = 660 * (variantList[i].pos - zoomSize[0]) / (zoomSize[1] - zoomSize[0]);
				var posY = 0;
				if (displayCoverage) {
					posY = 90 * (1 - (variantList[i].cov / maxCov));
					var lowerY = 90 * (1 - (covThreshold[0] / maxCov));
					var upperY = 90 * (1 - (covThreshold[1] / maxCov));
					nx.beginPath();
					nx.moveTo(0, lowerY);
					nx.lineTo(20, lowerY);
					nx.moveTo(0, upperY);
					nx.lineTo(20, upperY);
					nx.stroke();
				}
				else if (displayFrequency) {
					posY = 90 * (1 - (variantList[i].freq / 100));
					var lowerY = 90 * (1 - (freqThreshold[0] / 100));
					var upperY = 90 * (1 - (freqThreshold[1] / 100));
					nx.beginPath();
					nx.moveTo(0, lowerY);
					nx.lineTo(20, lowerY);
					nx.moveTo(0, upperY);
					nx.lineTo(20, upperY);
					nx.stroke();
				}
				else if (displayQuality) {
					posY = 90 * (1 - (variantList[i].qual / maxQual));
					var lowerY = 90 * (1 - (qualThreshold[0] / maxQual));
					var upperY = 90 * (1 - (qualThreshold[1] / maxQual));
					nx.beginPath();
					nx.moveTo(0, lowerY);
					nx.lineTo(20, lowerY);
					nx.moveTo(0, upperY);
					nx.lineTo(20, upperY);
					nx.stroke();
				}
				else if (displayBias) {
					posY = 90 * (1 - (variantList[i].strandBias));
					var lowerY = 90 * (1 - biasThreshold[0]);
					var upperY = 90 * (1 - biasThreshold[1]);
					nx.beginPath();
					nx.moveTo(0, lowerY);
					nx.lineTo(20, lowerY);
					nx.moveTo(0, upperY);
					nx.lineTo(20, upperY);
					nx.stroke();
				}
				var alt = variantList[i].varN;
				var ref = variantList[i].refN;
				var yCoor;
				if (ref == "G") {
					yCoor = 100;
				}
				else if (ref == "C") {
					yCoor = 115;
				}
				else if (ref == "T") {
					yCoor = 130;
				}
				else {
					yCoor = 145;
				}
				if (alt == "G") {
					ctx.strokeStyle = "#FF0000";
					ctx.fillStyle = "#FF0000";
				}
				else if (alt == "C") {
					ctx.strokeStyle = "#AFFF00";
					ctx.fillStyle = "#AFFF00";
				}
				else if (alt == "T") {
					ctx.strokeStyle = "#0000FF";
					ctx.fillStyle = "#0000FF";
				}
				else {
					ctx.strokeStyle = "00FFFF";
					ctx.fillStyle = "00FFFF";
				}
				ctx.fillText(alt, posX, yCoor);
				ctx.beginPath();
				ctx.moveTo(posX, posY);
				ctx.lineTo(posX, 90);
				ctx.stroke();
				cdx.beginPath();
				cdx.moveTo(posX, (100 - variantList[i].cadd)/100 * 25);
				cdx.lineTo(posX, 25);
				cdx.stroke();
				px.stroke();
				px.beginPath();
				px.moveTo(posX, (1 - variantList[i].phred) * 25);
				px.lineTo(posX, 25);
				px.stroke();
				y ++;
			}
			x ++;
		}
	}

	document.getElementById('filterRatio').value = y + "/" + x;
	displayCoverage = false;
	displayFrequency = false;
	displayQuality = false;
	displayBias = false;

}

function mousePressed(evt) {

	dragOn = true;
	var canvas = document.getElementById('mainTrack');
	var ctx = canvas.getContext('2d');
	ctx.strokeStyle = "#000000";
	ctx.fillStyle = "#000000";
	if (evt.clientX > Math.abs(evt.clientX - 660)) {
		dragSide = 1;
		ctx.fillRect(evt.clientX, 0, 660, 150);
	}
	else {
		dragSide = 0;
		ctx.fillRect(0, 0, evt.clientX, 150);
	}

}

function mouseDragged(evt) {

	if (dragOn) {
		var canvas = document.getElementById('mainTrack');
		var ctx = canvas.getContext('2d');
		ctx.strokeStyle = "#000000";
		ctx.fillStyle = "#000000";
		if (dragSide == 1) {
			ctx.fillRect(evt.clientX, 0, 660, 150);
		}
		else {
			ctx.fillRect(0, 0, evt.clientX, 150);
		}
	}

}

function mouseReleased(evt) {

	dragOn = false;
	var mouseX = zoomSize[0] + Math.round((evt.clientX / 660) * (zoomSize[1] - zoomSize[0]));
	zoomSize[dragSide] = mouseX;
	document.getElementById('start').value = zoomSize[0];
	document.getElementById('end').value = zoomSize[1];
	var c = document.getElementById('chrNum').value;
	document.getElementById('ucsc').style.visibility = "hidden";
	document.getElementById('ucsc').src="https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position=" + c + ":" + zoomSize[0] + "-" + zoomSize[1];

}

function mouseLeft(evt) {

	dragOn = false;
	draw();

}

function updateText(evt) {

	var n = evt.target.id;
	if (n == 'covRangeLower') {
		document.getElementById('covTextLower').value = evt.target.value;
		covThreshold[0] = evt.target.value;
		displayCoverage = true;
	}
	else if (n == 'covRangeUpper') {
		document.getElementById('covTextUpper').value = evt.target.value;
		covThreshold[1] = evt.target.value;
		displayCoverage = true;
	}
	else if (n == 'freqRangeLower') {
		document.getElementById('freqTextLower').value = evt.target.value;
		freqThreshold[0] = evt.target.value;
		displayFrequency = true;
	}
	else if (n == 'freqRangeUpper') {
		document.getElementById('freqTextUpper').value = evt.target.value;
		freqThreshold[1] = evt.target.value;
		displayFrequency = true;
	}
	else if (n == 'qualRangeLower') {
		document.getElementById('qualTextLower').value = evt.target.value;
		qualThreshold[0] = evt.target.value;
		displayQuality = true;
	}
	else if (n == 'qualRangeUpper') {
		document.getElementById('qualTextUpper').value = evt.target.value;
		qualThreshold[1] = evt.target.value;
		displayQuality = true;
	}
	else if (n == 'biasRangeLower') {
		document.getElementById('biasTextLower').value = evt.target.value;
		biasThreshold[0] = evt.target.value;
		displayBias = true;
	}
	else if (n == 'biasRangeUpper') {
		document.getElementById('biasTextUpper').value = evt.target.value;
		biasThreshold[1] = evt.target.value;
		displayBias = true;
	}
	draw();

}

function updateRange(evt) {

	var n = evt.target.id;
	if (n == 'covTextLower') {
		document.getElementById('covRangeLower').value = evt.target.value;
		covThreshold[0] = evt.target.value;
		displayCoverage = true;
	}
	else if (n == 'covTextUpper') {
		document.getElementById('covRangeUpper').value = evt.target.value;
		covThreshold[1] = evt.target.value;
		displayCoverage = true;
	}
	else if (n == 'freqTextLower') {
		document.getElementById('freqRangeLower').value = evt.target.value;
		freqThreshold[0] = evt.target.value;
		displayFrequency = true;
	}
	else if (n == 'freqTextUpper') {
		document.getElementById('freqRangeUpper').value = evt.target.value;
		freqThreshold[1] = evt.target.value;
		displayFrequency = true;
	}
	else if (n == 'qualTextLower') {
		document.getElementById('qualRangeLower').value = evt.target.value;
		qualThreshold[0] = evt.target.value;
		displayQuality = true;
	}
	else if (n == 'qualTextUpper') {
		document.getElementById('qualRangeUpper').value = evt.target.value;
		qualThreshold[1] = evt.target.value;
		displayQuality = true;
	}
	else if (n == 'biasTextLower') {
		document.getElementById('biasRangeLower').value = evt.target.value;
		biasThreshold[0] = evt.target.value;
		displayBias = true;
	}
	else if (n == 'biasTextUpper') {
		document.getElementById('biasRangeUpper').value = evt.target.value;
		biasThreshold[1] = evt.target.value;
		displayBias = true;
	}
	draw();

}

function changeGene(evt) {

	geneMode = geneNames.indexOf(evt.target.value);
	document.getElementById('ucsc').style.visibility = "hidden";
	document.getElementById('ucsc').src="https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position=chr" + geneChrs[geneMode] + ":" + geneEndPoints[geneMode][0] + "-" + geneEndPoints[geneMode][1];
	document.getElementById('chrNum').value = "chr" + geneChrs[geneMode];
	document.getElementById('start').value = geneEndPoints[geneMode][0];
	document.getElementById('end').value = geneEndPoints[geneMode][1];
	windowSize[0] = geneEndPoints[geneMode][0];
	windowSize[1] = geneEndPoints[geneMode][1];
	zoomSize[0] = geneEndPoints[geneMode][0];
	zoomSize[1] = geneEndPoints[geneMode][1];
}

function revertWindowStart(evt) {
	zoomSize[0] = windowSize[0];
	document.getElementById('ucsc').style.visibility = "hidden";
	document.getElementById('start').value = zoomSize[0];
	var c = document.getElementById('chrNum').value;
	document.getElementById('ucsc').src="https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position=" + c + ":" + zoomSize[0] + "-" + zoomSize[1];

}

function revertWindowEnd(evt) {
	zoomSize[1] = windowSize[1];
	document.getElementById('ucsc').style.visibility = "hidden";
	document.getElementById('end').value = zoomSize[1];
	var c = document.getElementById('chrNum').value;
	document.getElementById('ucsc').src="https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg38&position=" + c + ":" + zoomSize[0] + "-" + zoomSize[1];

}

function reset(evt) {
	covThreshold[0] = 0;
	covThreshold[1] = maxCov;
	freqThreshold[0] = 0;
	freqThreshold[1] = 100;
	qualThreshold[0] = 0;
	qualThreshold[1] = maxQual;
	biasThreshold[0] = 0;
	biasThreshold[1] = 1;
	document.getElementById('covRangeLower').value = 0;
	document.getElementById('covRangeUpper').value = maxCov;
	document.getElementById('covTextLower').value = 0;
	document.getElementById('covTextUpper').value = maxCov;
	document.getElementById('freqRangeLower').value = 0;
	document.getElementById('freqRangeUpper').value = 100;
	document.getElementById('freqTextLower').value = 0;
	document.getElementById('freqTextUpper').value = 100;
	document.getElementById('qualRangeLower').value = 0;
	document.getElementById('qualRangeUpper').value = maxQual;
	document.getElementById('qualTextLower').value = 0;
	document.getElementById('qualTextUpper').value = maxQual;
	document.getElementById('biasRangeLower').value = 0;
	document.getElementById('biasRangeUpper').value = 1;
	document.getElementById('biasTextLower').value = 0;
	document.getElementById('biasTextUpper').value = 1;
	draw();
}

document.getElementById('VCFUpload').addEventListener('change', readVCF, false);
document.getElementById('geneName').addEventListener('change', changeGene, false);

var doc = document.getElementById('ucsc');
doc.onload = function() {
		doc.style.visibility = "visible";
		draw();
}

document.getElementById('mainTrack').addEventListener('mousedown', mousePressed, false);
document.getElementById('mainTrack').addEventListener('mousemove', mouseDragged, false);
document.getElementById('mainTrack').addEventListener('mouseup', mouseReleased, false);
document.getElementById('mainTrack').addEventListener('mouseleave', mouseLeft, false);

document.getElementById('revertLeft').addEventListener('click', revertWindowStart, false);
document.getElementById('revertRight').addEventListener('click', revertWindowEnd, false);

document.getElementById('covRangeLower').addEventListener('input', updateText, false);
document.getElementById('covRangeUpper').addEventListener('input', updateText, false);
document.getElementById('covTextLower').addEventListener('change', updateRange, false);
document.getElementById('covTextUpper').addEventListener('change', updateRange, false);
document.getElementById('covDiv').addEventListener('mouseleave', draw, false);

document.getElementById('freqRangeLower').addEventListener('input', updateText, false);
document.getElementById('freqRangeUpper').addEventListener('input', updateText, false);
document.getElementById('freqTextLower').addEventListener('change', updateRange, false);
document.getElementById('freqTextUpper').addEventListener('change', updateRange, false);
document.getElementById('freqDiv').addEventListener('mouseleave', draw, false);

document.getElementById('qualRangeLower').addEventListener('input', updateText, false);
document.getElementById('qualRangeUpper').addEventListener('input', updateText, false);
document.getElementById('qualTextLower').addEventListener('change', updateRange, false);
document.getElementById('qualTextUpper').addEventListener('change', updateRange, false);
document.getElementById('qualDiv').addEventListener('mouseleave', draw, false);

document.getElementById('biasRangeLower').addEventListener('input', updateText, false);
document.getElementById('biasRangeUpper').addEventListener('input', updateText, false);
document.getElementById('biasTextLower').addEventListener('change', updateRange, false);
document.getElementById('biasTextUpper').addEventListener('change', updateRange, false);
document.getElementById('biasDiv').addEventListener('mouseleave', draw, false);

document.getElementById("resetFilters").addEventListener('click', reset, false);

// funtion memanggil button a href menu agar bisa mengakses main.js
$(function() {
	$('.navbar .navbar-nav a').not('.dropdown-toggle').bind('click', function(event) {
		var $anchor = $(this).attr('href'),
			$func   = $anchor.substr(1, $anchor.length);

		eval($func);
		return false;
	});
});

// buat variabel
var 
	cvs = document.getElementById("kanvas"), //variabel canvas
	ctx = cvs.getContext("2d"), //variabel dimensi, dimana 2d
	myImg = new Image(); //variabel gambar

// memanggil id imgFile pada modal upload
$(imgFile).change(function(){
	readURL(this);
});

//function read image
function readURL(input){
	if(input.files && input.files[0]){
		var reader =  new FileReader();
		reader.onload = function(e){
			$(myImg).attr("src", e.target.result);
		}
		reader.readAsDataURL(input.files[0]);
	}
}

// function untuk mengatur ukuran gambar saat ditampilkan di canvas
$(myImg).load(function(){
	$(kanvas).attr("width", myImg.width);
	$(kanvas).attr("height", myImg.height);
	// load original gambar
	ctx.drawImage(myImg, 0, 0);
	$("#inp-brightness").val(0);
	$("#inp-contrass").val(0);

	var imgData = ctx.getImageData(0, 0, cvs.width, cvs.height);
	$("#imgwidth").text("Lebar: "+imgData.width+"px");
	$("#imgheight").text("Tinggi: "+imgData.height+"px");
	$("#imglength").text("Panjang: "+imgData.data.length+"px");
	$(".imgtitle span").text("Gambar Asli");
});

//  function Reset gambar agar menjadi gambar asli
function imgReset(){
	$(kanvas).attr("width", myImg.width);
	$(kanvas).attr("height", myImg.height);
	$(".imgtitle span").text("Citra Asli");
	ctx.drawImage(myImg, 0, 0);
	$("#inp-brightness").val(0);
	$("#inp-contrass").val(0);
	$("#hist").hide();
}
// end function reset


// start funtion histogram
function histogram(){

	$("#hist").hide(); //untuk menutup histogram
	$("#bar-chart").css("width", "99.5%"); //lebar bar chart
	$("#bar-chart").css("height", "500px"); //tinggi bar chart

	// load original image
	// ctx.drawImage(myImg, 0, 0);

	// read img data
	var imgData = ctx.getImageData(0, 0, cvs.width, cvs.height);
	var hist = {"red":[], "green":[], "blue":[]};
	var pr = [], prrk = {"red":[], "green":[], "blue":[]};
	for(var n=0; n <= 255; n++){
		hist.red.push(0); hist.green.push(0); hist.blue.push(0);
		prrk.red.push(0); prrk.green.push(0); prrk.blue.push(0);
		pr.push(n/255);
	}

	// hitung nk
	for (var i = 0; i < imgData.data.length; i+=4) {
		for(var n=0; n <= 255; n++){
			if(imgData.data[i] == n) hist.red[n]++;
			if(imgData.data[i+1] == n) hist.green[n]++;
			if(imgData.data[i+2] == n) hist.blue[n]++;
		}

	}

	// Pr(rk)
	for(i=0; i<3; i++){
		for(var n=0; n <= 255; n++){
			if(i==0) prrk.red[n] = hist.red[n] / (imgData.data.length / 4);
			if(i==1) prrk.green[n] = hist.green[n] / (imgData.data.length / 4);
			if(i==2) prrk.blue[n] = hist.blue[n] / (imgData.data.length / 4);
		}
	}

	/* Bar Chart starts */

	var d1 = [];
	for (i = 0; i <= 255; i += 1)
		d1.push([i, hist.red[i]]);

	var d2 = [];
	for (i = 0; i <= 255; i += 1)
		d2.push([i, hist.green[i]]);

	var d3 = [];
	for (i = 0; i <= 255; i += 1)
		d3.push([i, hist.blue[i]]);

	var stack = 0, bars = true, lines = false, steps = false;

	function plotWithOptions() {
		$.plot($("#bar-chart"), [ d1, d2, d3 ], {
			series: {
				stack: stack,
				lines: { show: lines, fill: true, steps: steps },
				bars: { show: bars, barWidth: 0.8 }
			},
			grid: {
				borderWidth: 0, hoverable: true, color: "#777"
			},
			colors: ["#FF0000", "#00FF00", "#0000FF"],
			bars: {
				  show: true,
				  lineWidth: 0,
				  fill: true,
				  fillColor: { colors: [ { opacity: 0.9 }, { opacity: 0.8 } ] }
			}
		});
	}

	plotWithOptions();

	$(".stackControls input").click(function (e) {
		e.preventDefault();
		stack = $(this).val() == "With stacking" ? true : null;
		plotWithOptions();
	});
	$(".graphControls input").click(function (e) {
		e.preventDefault();
		bars = $(this).val().indexOf("Bars") != -1;
		lines = $(this).val().indexOf("Lines") != -1;
		steps = $(this).val().indexOf("steps") != -1;
		plotWithOptions();
	});

	/* Bar chart ends */

	$("#hist").slideDown('400');

}
// end function histogram

// start function Grayscale 
function imgGrayscale(){
	// load original gambar
	ctx.drawImage(myImg, 0, 0);

	// mengambil data gambar
	var imgData = ctx.getImageData(0, 0, cvs.width, cvs.height);

	// manipulation data gamabar
	for(var i=0; i < imgData.data.length; i +=4){
		// rumus yg digunakan menggunakan metode average
		// yaitu dengan menjumlahkan seluruh nilai R G B, kemudian dibagi 3, sehingga diperoleh nilai rata-rata dari R G B, nilai rata-rata itulah yang dapat dikatakan sebagai grayclase. 
		//rumus : Grayscale = (R + G + B) / 3 
		var gr = (imgData.data[i] + imgData.data[i+1] + imgData.data[i+2]) / 3;
		if(gr < 0) gr = 0;
		if(gr > 255) gr = 255;
		imgData.data[i] = gr; //R
		imgData.data[i+1] = gr; //G
		imgData.data[i+2] = gr; //B
	}

	// show manipulation
	ctx.putImageData(imgData, 0, 0);

	// set title
	$(".imgtitle span").text("Grayscale");
}
// end function Grayscale


// start funtion konvolusi
// konvolusi 3x3
function konvolusi3x3(inData, outData, width, height, kernel, alpha, invert, mono) {
	var idx, r, g, b, a,
		pyc, pyp, pyn,
		pxc, pxp, pxn,
		x, y,

		n = width * height * 4,

		k00 = kernel[0][0], k01 = kernel[0][1], k02 = kernel[0][2],
		k10 = kernel[1][0], k11 = kernel[1][1], k12 = kernel[1][2],
		k20 = kernel[2][0], k21 = kernel[2][1], k22 = kernel[2][2],

		p00, p01, p02,
		p10, p11, p12,
		p20, p21, p22;

	for (y=0;y<height;++y) {
		pyc = y * width * 4;
		pyp = pyc - width * 4;
		pyn = pyc + width * 4;

		if (y < 1) pyp = pyc;
		if (y >= width-1) pyn = pyc;

		for (x=0;x<width;++x) {
			idx = (y * width + x) * 4;

			pxc = x * 4;
			pxp = pxc - 4;
			pxn = pxc + 4;

			if (x < 1) pxp = pxc;
			if (x >= width-1) pxn = pxc;

			p00 = pyp + pxp;    p01 = pyp + pxc;    p02 = pyp + pxn;
			p10 = pyc + pxp;    p11 = pyc + pxc;    p12 = pyc + pxn;
			p20 = pyn + pxp;    p21 = pyn + pxc;    p22 = pyn + pxn;

			r = inData[p00] * k00 + inData[p01] * k01 + inData[p02] * k02
			  + inData[p10] * k10 + inData[p11] * k11 + inData[p12] * k12
			  + inData[p20] * k20 + inData[p21] * k21 + inData[p22] * k22;

			g = inData[p00 + 1] * k00 + inData[p01 + 1] * k01 + inData[p02 + 1] * k02
			  + inData[p10 + 1] * k10 + inData[p11 + 1] * k11 + inData[p12 + 1] * k12
			  + inData[p20 + 1] * k20 + inData[p21 + 1] * k21 + inData[p22 + 1] * k22;

			b = inData[p00 + 2] * k00 + inData[p01 + 2] * k01 + inData[p02 + 2] * k02
			  + inData[p10 + 2] * k10 + inData[p11 + 2] * k11 + inData[p12 + 2] * k12
			  + inData[p20 + 2] * k20 + inData[p21 + 2] * k21 + inData[p22 + 2] * k22;

			if (alpha) {
				a = inData[p00 + 3] * k00 + inData[p01 + 3] * k01 + inData[p02 + 3] * k02
				  + inData[p10 + 3] * k10 + inData[p11 + 3] * k11 + inData[p12 + 3] * k12
				  + inData[p20 + 3] * k20 + inData[p21 + 3] * k21 + inData[p22 + 3] * k22;
			} else {
				a = inData[idx+3];
			}

			if (mono) {
				r = g = b = (r + g + b) / 3;
			}

			if (invert) {
				r = 255 - r;
				g = 255 - g;
				b = 255 - b;
			}

			outData[idx] = r;
			outData[idx+1] = g;
			outData[idx+2] = b;
			outData[idx+3] = a;
		}
	}
}
// konvolusi 5x5
function konvolusi5x5(inData, outData, width, height, kernel, alpha, invert, mono) {
	var idx, r, g, b, a,
		pyc, pyp, pyn, pypp, pynn,
		pxc, pxp, pxn, pxpp, pxnn,
		x, y,

		n = width * height * 4,

		k00 = kernel[0][0], k01 = kernel[0][1], k02 = kernel[0][2], k03 = kernel[0][3], k04 = kernel[0][4],
		k10 = kernel[1][0], k11 = kernel[1][1], k12 = kernel[1][2], k13 = kernel[1][3], k14 = kernel[1][4],
		k20 = kernel[2][0], k21 = kernel[2][1], k22 = kernel[2][2], k23 = kernel[2][3], k24 = kernel[2][4],
		k30 = kernel[3][0], k31 = kernel[3][1], k32 = kernel[3][2], k33 = kernel[3][3], k34 = kernel[3][4],
		k40 = kernel[4][0], k41 = kernel[4][1], k42 = kernel[4][2], k43 = kernel[4][3], k44 = kernel[4][4],

		p00, p01, p02, p03, p04,
		p10, p11, p12, p13, p14,
		p20, p21, p22, p23, p24,
		p30, p31, p32, p33, p34,
		p40, p41, p42, p43, p44;

	for (y=0;y<height;++y) {
		pyc = y * width * 4;
		pyp = pyc - width * 4;
		pypp = pyc - width * 4 * 2;
		pyn = pyc + width * 4;
		pynn = pyc + width * 4 * 2;

		if (y < 1) pyp = pyc;
		if (y >= width-1) pyn = pyc;
		if (y < 2) pypp = pyp;
		if (y >= width-2) pynn = pyn;

		for (x=0;x<width;++x) {
			idx = (y * width + x) * 4;

			pxc = x * 4;
			pxp = pxc - 4;
			pxn = pxc + 4;
			pxpp = pxc - 8;
			pxnn = pxc + 8;

			if (x < 1) pxp = pxc;
			if (x >= width-1) pxn = pxc;
			if (x < 2) pxpp = pxp;
			if (x >= width-2) pxnn = pxn;

			p00 = pypp + pxpp;    p01 = pypp + pxp;    p02 = pypp + pxc;    p03 = pypp + pxn;    p04 = pypp + pxnn;
			p10 = pyp  + pxpp;    p11 = pyp  + pxp;    p12 = pyp  + pxc;    p13 = pyp  + pxn;    p14 = pyp  + pxnn;
			p20 = pyc  + pxpp;    p21 = pyc  + pxp;    p22 = pyc  + pxc;    p23 = pyc  + pxn;    p24 = pyc  + pxnn;
			p30 = pyn  + pxpp;    p31 = pyn  + pxp;    p32 = pyn  + pxc;    p33 = pyn  + pxn;    p34 = pyn  + pxnn;
			p40 = pynn + pxpp;    p41 = pynn + pxp;    p42 = pynn + pxc;    p43 = pynn + pxn;    p44 = pynn + pxnn;

			r = inData[p00] * k00 + inData[p01] * k01 + inData[p02] * k02 + inData[p03] * k04 + inData[p02] * k04
			  + inData[p10] * k10 + inData[p11] * k11 + inData[p12] * k12 + inData[p13] * k14 + inData[p12] * k14
			  + inData[p20] * k20 + inData[p21] * k21 + inData[p22] * k22 + inData[p23] * k24 + inData[p22] * k24
			  + inData[p30] * k30 + inData[p31] * k31 + inData[p32] * k32 + inData[p33] * k34 + inData[p32] * k34
			  + inData[p40] * k40 + inData[p41] * k41 + inData[p42] * k42 + inData[p43] * k44 + inData[p42] * k44;

			g = inData[p00+1] * k00 + inData[p01+1] * k01 + inData[p02+1] * k02 + inData[p03+1] * k04 + inData[p02+1] * k04
			  + inData[p10+1] * k10 + inData[p11+1] * k11 + inData[p12+1] * k12 + inData[p13+1] * k14 + inData[p12+1] * k14
			  + inData[p20+1] * k20 + inData[p21+1] * k21 + inData[p22+1] * k22 + inData[p23+1] * k24 + inData[p22+1] * k24
			  + inData[p30+1] * k30 + inData[p31+1] * k31 + inData[p32+1] * k32 + inData[p33+1] * k34 + inData[p32+1] * k34
			  + inData[p40+1] * k40 + inData[p41+1] * k41 + inData[p42+1] * k42 + inData[p43+1] * k44 + inData[p42+1] * k44;

			b = inData[p00+2] * k00 + inData[p01+2] * k01 + inData[p02+2] * k02 + inData[p03+2] * k04 + inData[p02+2] * k04
			  + inData[p10+2] * k10 + inData[p11+2] * k11 + inData[p12+2] * k12 + inData[p13+2] * k14 + inData[p12+2] * k14
			  + inData[p20+2] * k20 + inData[p21+2] * k21 + inData[p22+2] * k22 + inData[p23+2] * k24 + inData[p22+2] * k24
			  + inData[p30+2] * k30 + inData[p31+2] * k31 + inData[p32+2] * k32 + inData[p33+2] * k34 + inData[p32+2] * k34
			  + inData[p40+2] * k40 + inData[p41+2] * k41 + inData[p42+2] * k42 + inData[p43+2] * k44 + inData[p42+2] * k44;

			if (alpha) {
				a = inData[p00+3] * k00 + inData[p01+3] * k01 + inData[p02+3] * k02 + inData[p03+3] * k04 + inData[p02+3] * k04
				  + inData[p10+3] * k10 + inData[p11+3] * k11 + inData[p12+3] * k12 + inData[p13+3] * k14 + inData[p12+3] * k14
				  + inData[p20+3] * k20 + inData[p21+3] * k21 + inData[p22+3] * k22 + inData[p23+3] * k24 + inData[p22+3] * k24
				  + inData[p30+3] * k30 + inData[p31+3] * k31 + inData[p32+3] * k32 + inData[p33+3] * k34 + inData[p32+3] * k34
				  + inData[p40+3] * k40 + inData[p41+3] * k41 + inData[p42+3] * k42 + inData[p43+3] * k44 + inData[p42+3] * k44;
			} else {
				a = inData[idx+3];
			}

			if (mono) {
				r = g = b = (r + g + b) / 3;
			}

			if (invert) {
				r = 255 - r;
				g = 255 - g;
				b = 255 - b;
			}

			outData[idx] = r;
			outData[idx+1] = g;
			outData[idx+2] = b;
			outData[idx+3] = a;
		}
	}
}
// end function konvolusi

// start bagian output konvolusi
// start filte sharpen
function filter_sharpen(){

	var x, sx, sy, r, g, b, a, dstOff, srcOff, wt, cx, cy, scy, scx,
        weights = [0, -1, 0, -1, 5, -1, 0, -1, 0],
        katet = Math.round(Math.sqrt(weights.length)),
        half = (katet * 0.5) | 0,
        dstData = ctx.createImageData(cvs.width, cvs.height),
        dstBuff = dstData.data,
        srcBuff = ctx.getImageData(0, 0, cvs.width, cvs.height).data,
		y = cvs.height,
		mix = 0.9;
	while (y--) {
		x = cvs.width;
		while (x--) {
			sy = y;
			sx = x;
			dstOff = (y * cvs.width + x) * 4;
			r = 0;
			g = 0;
			b = 0;
			a = 0;

			for (cy = 0; cy < katet; cy++) {
				for (cx = 0; cx < katet; cx++) {
					scy = sy + cy - half;
					scx = sx + cx - half;

					if (scy >= 0 && scy < cvs.height && scx >= 0 && scx < cvs.width) {
						srcOff = (scy * cvs.width + scx) * 4;
						wt = weights[cy * katet + cx];

						r += srcBuff[srcOff] * wt;
						g += srcBuff[srcOff + 1] * wt;
						b += srcBuff[srcOff + 2] * wt;
						a += srcBuff[srcOff + 3] * wt;
					}
				}
			}

			dstBuff[dstOff] = r * mix + srcBuff[dstOff] * (1 - mix);
			dstBuff[dstOff + 1] = g * mix + srcBuff[dstOff + 1] * (1 - mix);
			dstBuff[dstOff + 2] = b * mix + srcBuff[dstOff + 2] * (1 - mix);
			dstBuff[dstOff + 3] = srcBuff[dstOff + 3];
		}
	}
	
	ctx.putImageData(dstData, 0, 0);
	

	// set title
	$(".imgtitle span").text("Sharpen");
}
// end filte sharpen
// start filte soften
function filter_soften(){

	// baca data gambar
	var imgData = ctx.getImageData(0, 0, cvs.width, cvs.height);
	var temp = imgData.data, iWidth = imgData.width, iHeight = imgData.height;

	var c = 1/9;
	konvolusi3x3(
		imgData.data, temp, iWidth, iHeight,
		[
		 [c, c, c],
		 [c, c, c],
		 [c, c, c]
		]
	);

	imgData.data = temp;

	// show manipulation
	ctx.putImageData(imgData, 0, 0);

	// set title
	$(".imgtitle span").text("Soften");
}
// end filte soften
// end bagian output konvolusi

// start threshold
function imgBinary(){

	// load original image
	ctx.drawImage(myImg, 0, 0);

	// get image data
	var imgData = ctx.getImageData(0, 0, cvs.width, cvs.height);

	// manipulation
	for(var i=0; i < imgData.data.length; i +=4){
		var gr = (imgData.data[i] + imgData.data[i+1] + imgData.data[i+2]) / 3;
		if(gr < 128) gr = 0;
		if(gr >= 128) gr = 255;
		imgData.data[i] = gr;
		imgData.data[i+1] = gr;
		imgData.data[i+2] = gr;
	}

	// show manipulation
	ctx.putImageData(imgData, 0, 0);

	// set title
	$(".imgtitle span").text("Gambar Threshold");

}
// end threshold

// start dilasi dan erosi
// start dilasi
function tepi_robert(){

	imgGrayscale();

	// read img data
	var imgData = ctx.getImageData(0, 0, cvs.width, cvs.height);
	var p, i, j, tempR = temp = imgData.data;
		iWidth 	= imgData.width, iHeight = imgData.height;

	// manipulation
	for(i = 0; i < (iHeight-1); i++) {
	  for(j = 0; j < (iWidth-1); j++) {

		  for (p=0; p<3; p++) { // p adalah penentu pixel r, g b

			  tempR[((iWidth * i) + j) * 4 + p] =
				  Math.abs( temp[((iWidth * i) + j) * 4 + p] - temp[((iWidth * (i+1)) + (j+1)) * 4 + p] ) +
				  Math.abs( temp[((iWidth * i) + (j+1)) * 4 + p] - temp[((iWidth * (i+1)) + j) * 4 + p] );

		  }

	  }
	}


	imgData.data = tempR;

	// show manipulation
	ctx.putImageData(imgData, 0, 0);

	// set title
	$(".imgtitle span").text("Deteksi Tepi Rebert");

}
// end dilasi

// start erosi
function tepi_prewitt(){

	imgGrayscale();

	// read img data
	var imgData = ctx.getImageData(0, 0, cvs.width, cvs.height);
	var p, i, j, tempR = temp = imgData.data;
		iWidth 	= imgData.width, iHeight = imgData.height;

	// manipulation
	for(i = 0; i < (iHeight-1); i++) {
	  for(j = 0; j < (iWidth-1); j++) {

		  for (p=0; p<3; p++) { // p adalah penentu pixel r, g b

			  var sx =
				  (temp[((iWidth * (i-1)) + (j-1)) * 4 + p] * (-1)) +
				  (temp[((iWidth * (i+0)) + (j-1)) * 4 + p] * (-1)) +
				  (temp[((iWidth * (i+1)) + (j-1)) * 4 + p] * (-1)) +

				  (temp[((iWidth * (i-1)) + (j+1)) * 4 + p] * (1)) +
				  (temp[((iWidth * (i+0)) + (j+1)) * 4 + p] * (1)) +
				  (temp[((iWidth * (i+1)) + (j+1)) * 4 + p] * (1));

			  var sy =
				  (temp[((iWidth * (i-1)) + (j-1)) * 4 + p] * (1)) +
				  (temp[((iWidth * (i-1)) + (j+0)) * 4 + p] * (1)) +
				  (temp[((iWidth * (i-1)) + (j+1)) * 4 + p] * (1)) +

				  (temp[((iWidth * (i+1)) + (j-1)) * 4 + p] * (-1)) +
				  (temp[((iWidth * (i+1)) + (j+0)) * 4 + p] * (-1)) +
				  (temp[((iWidth * (i+1)) + (j+1)) * 4 + p] * (-1));

			  var hit = 255 - (Math.abs(sx) + Math.abs(sy));
			  tempR[((iWidth * i) + j) * 4 + p] = hit;

		  }

	  }
	}


	imgData.data = tempR;

	// show manipulation
	ctx.putImageData(imgData, 0, 0);

	// set title
	$(".imgtitle span").text("Deteksi Tepi Prewitt");

}
// end erosi
// end dilasi dan erosi


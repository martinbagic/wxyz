let DATA; // gensurv, genrepr, phesurv, pherepr, deaths

let lifespan,
  bitsperlocus,
  maturationage,
  refresher,
  survivorship = [],
  svg = d3.select("svg"),
  rectSize = 12,
  colorant = d3.scaleSequential(d3.interpolateSpectral).domain([0, 1]),
  blockGap = 10,
  whiteGap = 1,
  maturationGap = 10,
  colorant2 = d3
    .scaleLinear()
    .domain([0, 1])
    .interpolate(d3.interpolateCubehelixLong)
    .range([d3.rgb("#0043de"), d3.rgb("#fc5a03")]), // fc5a03
  // colorant2 = d3.scaleSequential(d3.interpolateRdBu).domain([0, 1]),
  frames,
  barLength = 300,
  KMLength = 200,
  recordNum, // record to display
  deathColors = ["#36A700", "#EF6000", "#0044DE"]; // eco, gen, max
// deathColors = ["#700320", "#073467", "rgb(100,100,100)"]; //

// let genomeStruc = {
//   gensurv: [0, lifespan * bitsperlocus],
//   genrepr: [lifespan * bitsperlocus, lifespan * 2 * bitsperlocus],
//   genmuta: [lifespan * 2 * bitsperlocus, (lifespan * 2 + 1) * bitsperlocus],
//   phesurv: [
//     (lifespan * 2 + 1) * bitsperlocus,
//     (lifespan * 2 + 1) * bitsperlocus + lifespan,
//   ],
//   pherepr: [
//     (lifespan * 2 + 1) * bitsperlocus + lifespan,
//     (lifespan * 2 + 1) * bitsperlocus + lifespan * 2,
//   ],
// };

function transpose(a) {
  return Object.keys(a[0]).map(function (c) {
    return a.map(function (r) {
      return r[c];
    });
  });
}

let stack = d3
  .stack()
  .keys(d3.range(3))
  .order(d3.stackOrderNone)
  .offset(d3.stackOffsetNone);

function getStack(n) {
  let d = transpose([DATA.deatheco[n], DATA.deathgen[n], DATA.deathend[n]]);
  return stack(d);
}

function load(filename) {
  d3.json(`csvs/${filename}`, d3.autoType).then(function (d) {
    DATA = {
      gensurv: transpose(d.gensurv),
      genrepr: transpose(d.genrepr),
      phesurv: transpose(d.phesurv),
      pherepr: transpose(d.pherepr),
      deatheco: d.death_eco,
      deathend: d.death_end,
      deathgen: d.death_gen,
    };
    console.log(DATA);

    lifespan = d.lifespan;
    bitsperlocus = d.bitsperlocus;
    maturationage = d.maturationage;

    recordNum = DATA.deatheco.length - 1; // default record to display is the last one

    let sliderArea = svg.append("g").attr("transform", `translate(50,25)`);
    let plotArea = svg.append("g").attr("transform", `translate(50,50)`);

    frames = {
      gensurv: plotArea.append("g"),
      phesurv: plotArea
        .append("g")
        .attr(
          "transform",
          `translate(${rectSize * bitsperlocus + blockGap} 0)`
        ),
      pherepr: plotArea
        .append("g")
        .attr(
          "transform",
          `translate(${rectSize * (bitsperlocus + 1) + blockGap * 2} 0)`
        ),
      genrepr: plotArea
        .append("g")
        .attr(
          "transform",
          `translate(${rectSize * (bitsperlocus + 2) + blockGap * 3} 0)`
        ),
      deaths: plotArea
        .append("g")
        .attr(
          "transform",
          `translate(${rectSize * (bitsperlocus * 2 + 2) + blockGap * 4} 0)`
        ),
      kaplanmeier: plotArea
        .append("g")
        .attr(
          "transform",
          `translate(${rectSize * (bitsperlocus * 2 + 2) + blockGap * 4} 0)`
        ),
    };

    plotHeat(frames.gensurv, DATA.gensurv, bitsperlocus, recordNum);
    plotHeat(frames.phesurv, DATA.phesurv, 1, recordNum);
    plotHeat(frames.pherepr, DATA.pherepr, 1, recordNum);
    plotHeat(frames.genrepr, DATA.genrepr, bitsperlocus, recordNum);

    let slider = new Slider(
      sliderArea,
      rectSize * (bitsperlocus * 2 + 2) + blockGap * 3,
      recordNum
    );

    plotBars(frames.deaths, getStack(recordNum));
    console.log(DATA);

    // let deaths = transpose([DATA.deatheco, DATA.deathend, DATA.deathgen].map(a => a.map(x => d3.sum(x)))).map(x => d3.sum(x))
    // console.log(deaths);
    // let deaths;
    let deaths = transpose([
      DATA.deathgen,
      DATA.deatheco,
      DATA.deathend,
    ]).map((a) => transpose(a).map((x) => d3.sum(x)));
    for (let x of deaths) {
      let total = d3.sum(x);
      let cumsum = Array.from(d3.cumsum(x));
      cumsum.unshift(1);
      survivorship.push(cumsum.map((_) => (total - _) / total));
    }
    // plotKM(frames.kaplanmeier, recordNum);

    svg
      .attr("height", plotArea.node().getBBox().height + 50)
      .attr("width", 1000);

    // d3.select("#download").on("click", function () {
    //   download(d3.select("svg"));
    // });
  });
}

class Slider {
  constructor(frame, length, value) {
    this.scale = d3
      .scaleQuantize()
      .domain([0, length])
      .range(d3.range(value + 1));
    this.maxvalue = value;

    this.value = value;

    let slider = this;

    this.knobSize = rectSize;

    this.xlimit = length - this.knobSize;

    this.backLine = frame
      .append("line")
      .attr("x2", length)
      .attr("stroke", "#ddd")
      .attr("stroke-width", this.knobSize);

    function updateValue() {
      slider.value += 1;
      if (slider.value > slider.maxvalue) slider.value = 0;
      slider.knob.attr("x", slider.scale.invertExtent(slider.value)[0]);
      updatePlots(slider.value);
    }
    this.playing = false;
    this.timer = null;
    this.playButton = frame
      .append("rect")
      .attr("x", this.scale.invertExtent(this.value)[1] + blockGap)
      .attr("y", -this.knobSize / 2)
      .attr("height", this.knobSize)
      .attr("width", this.knobSize)
      .attr("fill", "red")
      .attr("cursor", "pointer")
      .on("click", function () {
        if (slider.playing) {
          slider.playing = false;
          slider.timer.stop();
          d3.select(this).attr("fill", "green");
        } else {
          slider.playing = true;
          slider.timer = d3.timer(updateValue, 1000); //slider.updateValue
          d3.select(this).attr("fill", "red");
        }
      });

    this.knob = frame
      .append("rect")
      .attr("x", this.scale.invertExtent(this.value)[0])
      .attr("y", -this.knobSize / 2)
      .attr("width", this.scale.invertExtent(0)[1])
      .attr("height", this.knobSize)
      .attr("fill", "rgb(0,68,222,0.8)")
      .attr("cursor", "grab")
      .call(
        d3
          .drag()
          .on("start", function (event) {
            d3.select(this).attr("opacity", 0.9);
          })
          .on("drag", function (event) {
            d3.select(this).style("cursor", "grabbing");
            d3.select("body").style("cursor", "grabbing");
            let x = event.x;
            if (x > slider.xlimit) x = slider.xlimit;
            if (x < 0) x = 0;

            slider.value = slider.scale(x);

            let xpos = slider.scale.invertExtent(slider.value)[0];
            d3.select(this).attr("x", xpos);

            updatePlots(slider.value);
          })
          .on("end", function (event) {
            d3.select("body").style("cursor", "default");
            d3.select(this).attr("opacity", 1).style("cursor", "grab");
          })
      );
  }
}

function updatePlots(recordNum) {
  for (let f of ["gensurv", "genrepr", "phesurv", "pherepr"]) {
    frames[f].selectAll("rect").attr("fill", (d) => colorant2(d[recordNum]));
  }

  plotBars(frames.deaths, getStack(recordNum));
  // plotKM(frames.kaplanmeier, recordNum);
}

function plotHeat(frame, data, bitsperlocus, recordNum) {
  function mouseover(event, d) {
    d3.select(this).style("opacity", 0.7);
    d3.select(".tooltip")
      .style("opacity", 1)
      .html(Math.round(d[recordNum] * 1000) / 1000)
      .style("left", event.pageX + "px")
      .style("top", event.pageY - 28 + "px")
      .style("color", colorant2(d[recordNum]));

    d3.select(".locuser")
      .style("opacity", 1)
      .style("color", "white")
      .html(Math.floor(d3.select(this).attr("i") / bitsperlocus))
      .style("left", event.pageX + "px")
      .style("top", event.pageY - 48 + "px");
  }

  function mouseout() {
    d3.select(this).style("opacity", 1);
    d3.selectAll(".tooltip, .locuser").style("opacity", 0);
  }

  function getY(_, i) {
    return (
      Math.floor(i / bitsperlocus) * rectSize +
      (i / bitsperlocus >= maturationage) * maturationGap
    );
  }

  frame
    .selectAll("rect")
    .data(data)
    .join("rect")
    .attr("i", (_, i) => i)
    .attr("height", rectSize - whiteGap)
    .attr("width", rectSize - whiteGap)
    .attr("x", (_, i) => (i % bitsperlocus) * rectSize)
    .attr("cursor", "pointer")
    .attr("y", getY)
    .attr("fill", (d) => colorant2(d[recordNum]))
    .on("mouseover", mouseover)
    .on("mouseout", mouseout);
}

function plotBars(frame, stack) {
  function mouseover(event, d) {
    d3.select(this).style("opacity", 0.7);
    d3.select(".tooltip")
      .style("opacity", 1)
      .html(d[1])
      .style("left", event.pageX + "px")
      .style("top", event.pageY - 28 + "px")
      .style("color", "white");
  }

  function mouseout() {
    d3.select(this).style("opacity", 1);
    d3.select(".tooltip").style("opacity", 0);
  }

  let maxBar = d3.max(stack[stack.length - 1], (x) => x[1]);
  let scale = d3.scaleLinear().domain([0, maxBar]).range([0, barLength]);

  frame
    .selectAll("g")
    .data(stack)
    .join("g")
    .attr("transform", (_, i) => `translate(${whiteGap * i}, 0)`)
    .selectAll("rect")
    .data((d, i) =>
      d.map((x) => {
        x.j = i;
        return x;
      })
    )
    .join("rect")
    .attr("cursor", "pointer")
    .attr("height", rectSize - 1)
    .attr("width", (d) => scale(d[1] - d[0]))
    .attr("y", (_, i) => rectSize * i + (i >= maturationage) * maturationGap)
    .attr("x", (d) => scale(d[0]))
    .attr("fill", (d) => deathColors[d.j]) // "#36A700", "#EF6000", "#0044DE"
    .on("mouseover", mouseover)
    .on("mouseout", mouseout);

  function addLines() {
    frame
      .selectAll("line")
      .data(d3.range(Math.floor(barLength / rectSize) + 1))
      .join("line")
      .attr("x1", (d) => rectSize * d)
      .attr("x2", (d) => rectSize * d)
      .attr("y1", 0)
      .attr("y2", lifespan * rectSize + maturationGap)
      .attr("stroke", "white")
      .attr("stroke-width", whiteGap);
  }

  // addLines();
}

function plotKM(frame, i) {
  var line = d3
    .line()
    .x((d) => d * KMLength)
    .y((_, i) => i * rectSize)
    .curve(d3.curveStepAfter);

  // frame
  //   .append("path")
  //   .datum(data)
  //   .attr("fill", "none")
  //   .attr("stroke", "black")
  //   .attr("stroke-width", 1.5)
  //   .attr("d", line);

  frame.selectAll("path").remove();

  let path = frame
    .selectAll("path")
    .data([survivorship[i]])
    .enter()
    .append("path")
    .attr("fill", "none")
    .attr("stroke", "black")
    .attr("stroke-width", 1.5)
    .attr("d", (d) => line(d))
    .exit();

  // path.exit().remove();
}

d3.select("#refreshme").on("click", function () {
  let removem = svg.selectAll("g");
  console.log(d3.select("input#filename").node().value);
  load(d3.select("input#filename").node().value);
  removem.remove();
});

// load("vizport.json");

// refresher = setInterval(function () {
//   let removem = svg.selectAll("g");
//   load("vizport.json");
//   removem.remove();
// }, 1000 * 10);
// clearInterval(refresher)

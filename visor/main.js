let data, // gensurv, genrepr, phesurv, pherepr, deaths
  survivorship = [],
  frames,
  recordNum; // record to display

const konst = {
  // preset
  rectSize: 12,
  barLength: 300,
  KMLength: 200,
  blockGap: 10,
  whiteGap: 1,
  maturationGap: 10,
  svg: d3.select("svg#svg1"),
  deathColors: ["#36A700", "#EF6000", "#0044DE"], // eco, gen, max
  colors: {
    blue: "#0044DE",
    green: "#36A700",
    orange: "#EF6000",
  },
  colorant: d3
    .scaleLinear()
    .domain([0, 1])
    .interpolate(d3.interpolateCubehelixLong)
    .range([d3.rgb("#0043de"), d3.rgb("#fc5a03")]),
  stack: d3
    .stack()
    .keys(d3.range(3))
    .order(d3.stackOrderNone)
    .offset(d3.stackOffsetNone),
  transpose: (a) => Object.keys(a[0]).map((c) => a.map((r) => r[c])),
  getStack: (i) =>
    konst.stack(
      konst.transpose([data.deatheco[i], data.deathgen[i], data.deathend[i]])
    ),
  // from input
  bitsperlocus: null,
  maturationage: null,
  maxlifespan: null,
};

const plotter = {
  KaplanMeier: () => {
    let line = d3
      .line()
      .x((d) => d * konst.KMLength)
      .y((_, i) => i * konst.rectSize)
      .curve(d3.curveStepAfter);

    frames.kaplanmeier.selectAll("path").remove();

    // line
    frames.kaplanmeier
      .selectAll("path")
      .data([survivorship[recordNum]])
      .enter()
      .append("path")
      .attr("fill", "none")
      .attr("stroke", "black")
      .attr("stroke-width", 1.5)
      .attr("d", (d) => line(d))
      .exit();
  },
  HeatTiles: (typ, bitsperlocus = konst.bitsperlocus) => {
    let mouseover = (event, d) => {
      plotter.LocusTimeline(d3.select(event.target).attr("i"));
      d3.select(event.target).style("opacity", 0.7);
      d3.select(".tooltip")
        .style("opacity", 1)
        .html(Math.round(d[recordNum] * 1000) / 1000)
        .style("left", event.pageX + "px")
        .style("top", event.pageY - 28 + "px")
        .style("color", konst.colorant(d[recordNum]));

      d3.select(".locuser")
        .style("opacity", 1)
        .style("color", "white")
        .html(Math.floor(d3.select(event.target).attr("i") / bitsperlocus))
        .style("left", event.pageX + "px")
        .style("top", event.pageY - 48 + "px");
    };

    let mouseout = (event) => {
      d3.select(event.target).style("opacity", 1);
      d3.selectAll(".tooltip, .locuser").style("opacity", 0);

      d3.selectAll("div.locusTimelines > svg.delete").remove();
    };

    let getY = (_, i) =>
      Math.floor(i / bitsperlocus) * konst.rectSize +
      (i / bitsperlocus >= konst.maturationage) * konst.maturationGap;

    frames[typ]
      .selectAll("rect")
      .data(data[typ])
      .join(
        (enter) =>
          enter
            .append("rect")
            .attr("i", (_, i) => i)
            .attr("height", konst.rectSize - konst.whiteGap)
            .attr("width", konst.rectSize - konst.whiteGap)
            .attr("x", (_, i) => (i % bitsperlocus) * konst.rectSize)
            .attr("cursor", "pointer")
            .attr("y", getY)
            .attr("fill", (d) => konst.colorant(d[recordNum]))
            .on("mouseover", mouseover)
            .on("mouseout", mouseout)
            .on("click", (event) => {
              plotter.LocusTimeline(d3.select(event.target).attr("i"), true);
            }),
        (update) => update.attr("fill", (d) => konst.colorant(d[recordNum]))
      );
  },
  DeathBars: () => {
    let stack = konst.getStack(recordNum);
    let mouseover = (event, d) => {
      d3.select(event.target).style("opacity", 0.7);
      d3.select(".tooltip")
        .style("opacity", 1)
        .html(d[1])
        .style("left", event.pageX + "px")
        .style("top", event.pageY - 28 + "px")
        .style("color", "white");
    };

    let mouseout = (event) => {
      d3.select(event.target).style("opacity", 1);
      d3.select(".tooltip").style("opacity", 0);
    };

    let maxBar = d3.max(stack[stack.length - 1], (x) => x[1]);
    let scale = d3
      .scaleLinear()
      .domain([0, maxBar])
      .range([0, konst.barLength]);

    frames.deaths
      .selectAll("g")
      .data(stack)
      .join("g")
      .attr("transform", (_, i) => `translate(${konst.whiteGap * i}, 0)`)
      .selectAll("rect")
      .data((d, i) =>
        d.map((x) => {
          x.j = i;
          return x;
        })
      )
      .join("rect")
      .attr("cursor", "pointer")
      .attr("height", konst.rectSize - 1)
      .attr("width", (d) => scale(d[1] - d[0]))
      .attr(
        "y",
        (_, i) =>
          konst.rectSize * i + (i >= konst.maturationage) * konst.maturationGap
      )
      .attr("x", (d) => scale(d[0]))
      .attr("fill", (d) => konst.deathColors[d.j]) // "#36A700", "#EF6000", "#0044DE"
      .on("mouseover", mouseover)
      .on("mouseout", mouseout);
  },
  LocusTimeline: (tilej, keep = false) => {
    // let deaths = konst
    //   .transpose([data.deathgen, data.deatheco, data.deathend])
    //   .map((a) => konst.transpose(a).map((x) => d3.sum(x)));
    // for (let x of deaths) {
    //   let total = d3.sum(x);
    //   let cumsum = Array.from(d3.cumsum(x));
    //   cumsum.unshift(1);
    //   survivorship.push(cumsum.map((_) => (total - _) / total));
    // }

    let plotH = 100;
    let plotW = 400;

    d3.select("div.locusTimelines").style("height", `${plotH}px`);
    d3.select("#locusTimelineWrapper").style("width", `${plotW}px`);

    let frameSvg = d3
      .select(keep ? "div.locusTimelinesSaved" : "div.locusTimelines")
      .append("svg")
      .classed("delete", keep ? false : true)
      .attr("width", plotW)
      .attr("height", plotH);

    let frame = frameSvg.append("g"); //.attr("transform", "translate(20, 5)");

    if (keep) {
      function save() {
        d3.selectAll("div.locusTimelinesSaved > svg").attr("id", null);
        frameSvg.attr("id", "tobesaved").attr("filename", "bazuka");
        SaveSVG.save("tobesaved", "./styles.css", "");
        return Promise.resolve(true);
      }

      function remove() {
        frameSvg.remove();
      }
      frame
        .append("rect")
        .attr("height", 20)
        .attr("width", 20)
        .on("click", remove);

      frame
        .append("rect")
        .attr("height", 20)
        .attr("width", 20)
        .attr("y", 30)
        .style("fill", "green")
        .on("click", save);

      frame
        .append("rect")
        .attr("height", 20)
        .attr("width", 20)
        .attr("y", 60)
        .style("fill", "pink")
        .on("click", () => {
          save().then(() => remove());
        });
    }

    let bitj = tilej % konst.bitsperlocus;
    let locusj = Math.floor(tilej / konst.bitsperlocus);

    let timelineData = data.gensurv.slice(
      locusj * konst.bitsperlocus,
      (locusj + 1) * konst.bitsperlocus
    );

    let wUnit = plotW / timelineData[0].length;

    let line = d3
      .line()
      .x((_, i) => wUnit * i)
      .y((d) => (1 - d) * plotH)
      .curve(d3.curveCardinal);
    // .curve(d3.curveStepAfter);

    // vertical line
    // frames.locusarea
    //   .append("line")
    //   .attr("y1", 0)
    //   .attr("y2", plotH)
    //   .attr("x1", recordNum * wUnit)
    //   .attr("x2", recordNum * wUnit)
    //   .style("stroke", "#333")
    //   .style("stroke-width", 2);

    // add axes
    // var xscale = d3
    //   .scaleLinear()
    //   .domain([0, timelineData[0].length])
    //   .range([0, plotW]);
    // var yscale = d3.scaleLinear().domain([1, 0]).range([0, plotH]);
    // var x_axis = d3.axisBottom().scale(xscale).ticks(5);
    // var y_axis = d3.axisLeft().scale(yscale).ticks(3);
    // frame
    //   .append("g")
    //   .attr("transform", `translate(${20} ${plotH})`)
    //   .call(x_axis);
    // frame.append("g").attr("transform", `translate(${20} ${0})`).call(y_axis);

    // lines
    frame
      .append("g")
      // .attr("transform", `translate(20 0)`)
      .selectAll("path")
      .data(timelineData)
      .enter()
      .append("path")
      .attr("fill", "none")
      .attr("stroke", (d, i) =>
        i == bitj ? konst.colorant(d[recordNum]) : "#333"
      )
      .attr("stroke-width", (_, i) => (i == bitj ? 2 : 0.5))
      .attr("stroke-dasharray", (_, i) => (i == bitj ? "0" : "3,6"))
      .style("stroke-linecap", "rounded")
      .attr("d", (d) => line(d))
      .exit();

    refitSVG();
  },
  replot: () => {
    plotter.HeatTiles("gensurv");
    plotter.HeatTiles("phesurv", 1);
    plotter.HeatTiles("pherepr", 1);
    plotter.HeatTiles("genrepr");
    plotter.DeathBars();
    plotter.KaplanMeier();

    // make sure svg height fits
    refitSVG();
  },
};

const refitSVG = () => {
  konst.svg
    .attr("height", frames.locusarea.node().getBoundingClientRect().bottom)
    .attr("width", frames.plotarea.node().getBoundingClientRect().right);
};

const makeSlider = (frame, length, value) => {
  // attrs
  let self = {
    maxvalue: value,
    value: value,
    knobSize: konst.rectSize,
    playing: false,
    timer: null,
  };

  self.xlimit = length - self.knobSize;
  self.scale = d3
    .scaleQuantize()
    .domain([0, self.xlimit])
    .range(d3.range(self.value + 1)); // max record index + 1 because d3.range stop value is exclusive

  // backline
  frame
    .append("line")
    .attr("x2", length)
    .attr("stroke", "#ddd")
    .attr("stroke-width", self.knobSize);

  // play button
  frame
    .append("rect")
    .attr("x", length + konst.blockGap)
    .attr("y", -self.knobSize / 2)
    .attr("height", self.knobSize)
    .attr("width", self.knobSize)
    .attr("fill", "gray")
    .attr("cursor", "pointer")
    .on("click", (event) => {
      if (self.playing) {
        self.playing = false;
        self.timer.stop();
        d3.select(event.target).attr("fill", "gray");
      } else {
        self.playing = true;
        self.timer = d3.interval(updateValue, 5000 / self.maxvalue);
        d3.select(event.target).attr("fill", "red");
      }
    });

  // current rec num
  let recText = frame
    .append("text")
    .attr("x", length + konst.blockGap * 2 + self.knobSize)
    .attr("y", 4)
    .style("font-size", 12)
    .style("font-family", "Monaco");

  let updateRecText = () => recText.text(`(0 - ${self.maxvalue}) ${recordNum}`);
  updateRecText();

  // knob
  let knob = frame
    .append("rect")
    .attr("x", self.scale.invertExtent(self.value)[1])
    .attr("y", -self.knobSize / 2)
    // .attr("width", scale.invertExtent(0)[1])
    .attr("width", self.knobSize)
    .attr("height", self.knobSize)
    .attr("fill", "rgb(0,68,222,0.8)")
    .attr("cursor", "grab")
    .call(
      d3
        .drag()
        .on("start", function (event) {
          console.log(recordNum);
          d3.select(this).attr("opacity", 0.9);
        })
        .on("drag", function (event) {
          // console.log(this);
          // d3.select(event.target).style("cursor", "grabbing");
          d3.select("body").style("cursor", "grabbing");
          let x = event.x;
          if (x > self.xlimit) x = self.xlimit;
          if (x < 0) x = 0;

          self.value = self.scale(x);

          // let xpos = self.scale.invertExtent(self.value)[0];
          d3.select(this).attr("x", x);

          recordNum = self.value;
          updateRecText();
          plotter.replot();
        })
        .on("end", function (event) {
          d3.select("body").style("cursor", "default");
          d3.select(this).attr("opacity", 1).style("cursor", "grab");
        })
    );

  let updateValue = () => {
    self.value += 1;
    if (self.value > self.maxvalue) self.value = 0;
    knob.attr("x", self.scale.invertExtent(self.value)[0]);
    recordNum = self.value;
    updateRecText();
    plotter.replot();
  };
};

const setVars = {
  setData: (raw) => {
    data = {
      gensurv: konst.transpose(raw.gensurv),
      genrepr: konst.transpose(raw.genrepr),
      phesurv: konst.transpose(raw.phesurv),
      pherepr: konst.transpose(raw.pherepr),
      deatheco: raw.death_eco,
      deathend: raw.death_end,
      deathgen: raw.death_gen,
    };
  },
  setRecordNum: () => {
    // default record to display is the last one
    recordNum = data.deatheco.length - 1;
  },

  setFrames: () => {
    let sliderArea = konst.svg
      .append("g")
      .attr("transform", `translate(50,25)`);
    let plotArea = konst.svg.append("g").attr("transform", "translate(50,50)");

    let locusArea = konst.svg
      .append("g")
      .attr(
        "transform",
        `translate(50,${
          konst.rectSize * konst.maxlifespan + konst.blockGap + 75
        })`
      )
      .append("svg")
      .attr("id", "locusTimeline");

    let getFrame = (locilens, rectlens, blockgaps) =>
      plotArea
        .append("g")
        .attr(
          "transform",
          `translate(${
            konst.rectSize * (konst.bitsperlocus * locilens + rectlens) +
            konst.blockGap * blockgaps
          } 0)`
        );

    frames = {
      gensurv: getFrame(0, 0, 0),
      phesurv: getFrame(1, 0, 1),
      genrepr: getFrame(1, 2, 3),
      pherepr: getFrame(1, 1, 2),
      deaths: getFrame(2, 2, 4),
      kaplanmeier: getFrame(2, 2, 4),
      sliderarea: sliderArea,
      plotarea: plotArea,
      locusarea: locusArea,
    };
  },

  setSurvivorship: () => {
    let deaths = konst
      .transpose([data.deathgen, data.deatheco, data.deathend])
      .map((a) => konst.transpose(a).map((x) => d3.sum(x)));
    for (let x of deaths) {
      let total = d3.sum(x);
      let cumsum = Array.from(d3.cumsum(x));
      cumsum.unshift(1);
      survivorship.push(cumsum.map((_) => (total - _) / total));
    }
  },
};

async function load(filename) {
  let rawData = await d3.json(`csvs/${filename}`, d3.autoType);
  konst.bitsperlocus = rawData.bitsperlocus;
  konst.maturationage = rawData.maturationage;
  konst.maxlifespan = rawData.lifespan;

  setVars.setData(rawData);
  setVars.setRecordNum();
  setVars.setFrames();
  setVars.setSurvivorship();

  makeSlider(
    frames.sliderarea,
    konst.rectSize * (konst.bitsperlocus * 2 + 2) + konst.blockGap * 3,
    recordNum // max recordNum
  );

  plotter.replot();
}

function refresh() {
  let filename = d3.select("input#filename").node().value;
  console.log(`Refreshing from ${filename}.`);
  konst.svg.selectAll("g").remove();
  load(filename);
}

// detect refresh key
{
  let keyIsDown = false;
  d3.select("body")
    .on("keydown", (event) => {
      if (keyIsDown) return;
      keyIsDown = true;
      if (event.key == "Enter") refresh();
    })
    .on("keyup", () => {
      keyIsDown = false;
    });
}
// load me button
d3.select("#refreshme").on("click", () => refresh());

// keep refreshing
{
  let keepRefreshing = null;
  d3.select("#keepRefreshing").on("click", function () {
    if (keepRefreshing == null) {
      console.log("Refreshing initiated.");
      keepRefreshing = setInterval(refresh, 1000 * 10);
      d3.select(this).classed("buttonOn", true);
    } else {
      console.log("Refreshing aborted.");
      clearInterval(keepRefreshing);
      keepRefreshing = null;
      d3.select(this).classed("buttonOn", false);
    }
  });
}

// download button
d3.select("#download").on("click", function () {
  console.log("Saving SVG.");
  SaveSVG.save("ana", "./styles.css", "");
});

// initial load
load("example.json");

// make deleteAll button delete all svgs

d3.select("#locusTimelineWrapper > .deleteAll").on("click", () => {
  d3.selectAll(".locusTimelinesSaved > svg").remove();
});

const transpose = (a) => Object.keys(a[0]).map((c) => a.map((r) => r[c]));

const conf = {
  svg: d3.select("svg#svg1"),
};

async function load(filename) {
  let d = await d3.json(`csvs/${filename}`, d3.autoType);

  DATA = {
    gensurv: transpose(d.gensurv),
    genrepr: transpose(d.genrepr),
    phesurv: transpose(d.phesurv),
    pherepr: transpose(d.pherepr),
    deatheco: d.death_eco,
    deathend: d.death_end,
    deathgen: d.death_gen,
  };
  //   console.log("DATA", DATA);
  console.log(d);
}

load("visor.json");

function plotAll() {
  function plotTiles(frame, data, stagei) {
    frame
      .selectAll("rect")
      .data(data)
      .join(
        (enter) => enter.append("rect").attr("height", 10).attr("width", 10).attr("x"),
        (update) => update.attr("fill", (d) => colorant2(d[stagei]))
      );
  }
}

let svg = d3.select("svg")


let data;

let vizportStatus = d3.select("#vizportStatus")

let heatmapGroup = svg.append("g").classed("heatmapGroup", true).attr("transform", "translate(0,0)")

let lineplotGroup = svg.append("g").classed("lineplotGroup", true)
let barplotGroup = svg.append("g").classed("barplotGroup", true)


function loadAndPlot(fileName) {
    d3.csv(`csvs/${fileName}`, d3.autoType).then(
        function (d) {
            data = d;
            vizportStatus.text("");
            d3.select("title").text("AEGIS VizPort: " + fileName);
            plot(data);
        },
        function () {
            vizportStatus.text("failed");
            setTimeout(_ => vizportStatus.text(""), 5000)
        }
    )
}


let vizportInput = d3.select("input#vizportFile")
    .on("keypress", function (event) {
        if (event.key === "Enter") {
            vizportStatus.text("loading");
            d3.selectAll("svg > .heatmapGroup > g").remove();
            d3.selectAll("svg > .lineplotGroup > *").remove();
            let fileName = d3.select("input#vizportFile").property("value");
            console.log(fileName)
            loadAndPlot(fileName)

        }
    })

vizportInput.node().focus();




function plot(data) {

    let arrays = {
        g: [],
        e: [],
        a: [],
        cod_: [],
    }

    let getLastNElements = 50;


    for (col of data.columns) {
        for (kind in arrays) {
            if (col.startsWith(kind)) {
                arrays[kind].push(data.map(_ => _[col]));
                break;
            }
        }
    };

    for (kind in arrays) {
        let len = arrays[kind][0].length - getLastNElements;
        if (len > 0) {
            arrays[kind] = arrays[kind].map(_ => _.slice(len))
        }
    }

    let bitsPerLocus = Math.floor(arrays["g"].length / arrays["e"].length);


    let get2DMax = array2d => d3.max(array2d, _ => d3.max(_))


    let scaleRange = ["#0043de", "#0090de", "#8dde00", "#e6db0b", "#fc5a03"]//,"#d94d03"];

    let getScale = (array2d, max = null, range = null) => d3.scaleLinear()
        .domain(d3.ticks(0, max == null ? get2DMax(array2d) : max, scaleRange.length))
        .range(range == null ? scaleRange : range)

    let xOffset = 0;
    let rwid = 6;
    let frame;

    let scales = {
        g: getScale(arrays["g"], 1),
        // e: getScale(arrays["e"], 1.25, ["#0043de","#0090de","#8dde00","#fc5a03", "#ffa97a"]),
        e: getScale(arrays["e"], 1),
        a: getScale(arrays["a"]),
        cod_: getScale(arrays["cod_"]),
    }

    for (arrayk in arrays) {
        console.log(arrayk);
        frame = heatmapGroup.append("g").attr("transform", `translate(${xOffset},0)`);
        heatmap(frame, arrays[arrayk], rwid, scales[arrayk], arrayk == "g" ? bitsPerLocus : 10, arrayk);
        xOffset += arrays[arrayk][0].length * rwid + 3;
    }


    let lineplotHeight = 100;
    let lineplotWidth = arrays["cod_"][0].length * 10;

    let svgHeight = rwid * 1.1 * d3.max(Object.keys(arrays).map(key => arrays[key].length));
    let svgWidth = arrays.g[0].length * rwid * Object.keys(arrays).length + 30;

    svg.attr("height", svgHeight)
    svg.attr("width", svgWidth + lineplotWidth + 100)

    lineplotGroup.attr("transform", `translate(${svgWidth},0)`)

    lineplot(lineplotGroup, arrays["cod_"], lineplotWidth, lineplotHeight);

    // svgWidth += 100;
    // barplotGroup.attr("transform", `translate(${svgWidth},0)`)
    // barplot(barplotGroup, arrays["cod_"], rwid);

}

loadAndPlot(d3.select("input#vizportFile").attr("value"))

let svg = d3.select("svg")


let data;

let vizportStatus = d3.select("#vizportStatus")
// let dataLoadingStatus = svg.append("text")
//     .attr("x", 20)
//     .attr("y", 20)


let heatmapGroup = svg.append("g").classed("heatmapGroup", true).attr("transform", "translate(0,0)")


// d3.select("input#vizportFileButton")
//     .on("click", function () {
//         // d3.selectAll(".heatmapGroup.g").delete();
//         dataLoadingStatus.text("loading");
//         let fileName = d3.select("input#vizportFile").property("value");
//         console.log(fileName)
//         d3.csv(fileName).then(
//             function (d) {
//                 data = d;
//                 dataLoadingStatus.text("done");
//                 // setTimeout(_ => dataLoadingStatus.text(""), 1500)
//                 plot(data);
//             },
//             function () {
//                 dataLoadingStatus.text("failed");
//                 setTimeout(_ => dataLoadingStatus.text(""), 5000)
//             }
//         )
//         // console.log(x)
//     })

let vizportInput = d3.select("input#vizportFile")
    .on("keypress", function (event) {
        if (event.key === "Enter") {
            vizportStatus.text("loading");
            console.log(d3.selectAll("svg > .heatmapGroup"));
            d3.selectAll("svg > .heatmapGroup > g").remove();
            let fileName = d3.select("input#vizportFile").property("value");
            console.log(fileName)
            d3.csv(fileName).then(
                function (d) {
                    data = d;
                    console.log(data)
                    vizportStatus.text("");
                    plot(data);
                },
                function () {
                    vizportStatus.text("failed");
                    setTimeout(_ => vizportStatus.text(""), 5000)
                }
            )
        }
    })

vizportInput.node().focus();



// let frames = {
//     ages: svg.append("g").attr("transform", "translate(0, 0)"),
//     heatma: svg.append("g").attr("transform", "translate(400, 400)"),
//     heatma: svg.append("g").attr("transform", "translate(400, 0)"),
// }


function plot(data) {

    let arrays = {
        g: [],
        e: [],
        a: [],
        // cod_: [],
    }

    for (col of data.columns) {
        for (kind in arrays) {
            if (col.startsWith(kind)) {
                arrays[kind].push(data.map(_ => _[col]));
                break;
            }
        }
    };

    let getScale = array2d => d3.scaleSequential()
        .domain([0, d3.max(array2d, _ => d3.max(_))])
        .interpolator(d3.scaleSequential(["blue", "red"]));

    let xOffset = 0;
    let rwid = 5;
    let frame;

    for (arrayk in arrays) {
        console.log(arrayk);
        frame = heatmapGroup.append("g").attr("transform", `translate(${xOffset},0)`);
        heatmap(frame, arrays[arrayk], rwid, getScale(arrays[arrayk]));
        xOffset += arrays[arrayk][0].length * rwid;
    }

    console.log(d3.max(Object.keys(arrays).map(key => arrays[key].length)))

    svg.attr("height", rwid * d3.max(Object.keys(arrays).map(key => arrays[key].length)))
    svg.attr("width", arrays.g[0].length * rwid * Object.keys(arrays).length)

}


// d3.csv(fileName).then(
//     function (d) {
//         data = d;
//         // dataLoadingStatus.text("done");
//         console.log(data);
//         // console.log(getcol(data, "birthdays"))
//         setTimeout(_ => dataLoadingStatus.text(""), 1500)
//         plot(data);
//     },
//     function () {
//         dataLoadingStatus.text("failed");
//         setTimeout(_ => dataLoadingStatus.text(""), 5500)
//     }
// )
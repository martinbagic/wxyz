const zip = (arr1, arr2) => arr1.map((k, i) => [k, arr2[i]]);

let tooltip = d3.select("div.tooltip")
let locuser = d3.select("div.locuser")

function heatmap(frame, data, rwid, scale, bitsperlocus, arrayk) {

    let dataSpan = d3.select("span#locus")

    let lines = frame
        .selectAll("g")
        .data(data)
        .join("g")
        .attr("transform", (_, i) => `translate(0, ${i * rwid + Math.floor(i / bitsperlocus) * 1})`)
        .attr("i", (_, i) => i)
        .attr("arrayk", arrayk)
        .on("mouseover", function (event, d) {
            let self = d3.select(this);
            self
                .style("opacity", 0.9)

            let diviz = arrayk == "g" ? bitsperlocus : 1;

            locuser.style("opacity", 1)
                .style("color", "white")
                .html(Math.floor(self.attr("i") / diviz))
                .style("left", (event.pageX) + "px")
                .style("top", (event.pageY - 48) + "px")
        })
        .on("mouseout", function (event, d) {
            d3.select(this)
                .style("opacity", 1)
            locuser
                .style("opacity", 0);
        });


    lines
        .selectAll("rect")
        .data(d => d)
        .join("rect")
        .attr("width", rwid - 0.3) // -0.5
        .attr("height", rwid - 0.3) // -0.5
        .attr("x", (_, i) => rwid * i)
        .attr("order", (_, i) => i)
        // .style("stroke-width", 0.1)
        // .style("stroke", "#000")
        .style("rx", 0.6)
        .style("ry", 0.6)
        .style("fill", d => scale(d))
        .on("mouseover", function (event, d, i) {
            let self = d3.select(this);
            self
                .style("opacity", 0.7)
            tooltip.transition()
                .duration(200)
                .style("opacity", 1);
            tooltip.html(Math.round(d * 1000) / 1000)
                .style("left", (event.pageX) + "px")
                .style("top", (event.pageY - 28) + "px")
                .style("color", scale(d));
        })
        .on("mouseout", function (event, d) {
            d3.select(this)
                .style("opacity", 1)
            tooltip.transition()
                .duration(500)
                .style("opacity", 0);
        });


}

function lineplot(frame, deaths, width, height) {

    let margin = 5

    let maxDeath = d3.max(deaths, _ => d3.max(_));


    let xScale = d3.scaleLinear()
        .domain([0, deaths[1].length - 1])
        .range([margin, width - margin]);

    let yScale = d3.scaleLinear()
        .domain([0, maxDeath])
        .range([height - margin, margin]);

    // frame.append("g").call(d3.axisLeft(yScale))
    // frame.append("g")
    //     .attr("transform", `translate(0, ${height})`)
    //     .call(d3.axisBottom(xScale))

    frame.append("rect")
        .style("height", height - margin / 2)
        .style("width", width - margin / 2)
        .style("stroke-width", 0.2)
        .style("stroke", "black")
        .style("fill", "transparent")
        .attr("x", margin / 4)
        .attr("y", margin / 4)


    let line = d3.line()
        .x(function (_, i) { return xScale(i); })
        .y(function (d) { return yScale(d); })
        .curve(d3.curveMonotoneX)

    // var dataset = d3.range(10).map(function (_) { return { "y": d3.randomUniform(1)() } })



    let lines = frame.selectAll("lines")
        .data(deaths)
        .enter()
        .append("g")


    let lineColors = d3.scaleSequential()
        .domain([0, deaths.length])
        .interpolator(d3.interpolateSpectral)

    lines.append("path")
        .attr("d", function (d) { return line(d); })
        // .style("stroke", "#0043de")
        .style("stroke-width", 1)
        .style("line-style", "dotted")
        .style("opacity", 0.8)
        .style("fill", "none")
        .style("stroke", (_, i) => lineColors(i))

    // let getScale = array2d => d3.scaleSequential()
    //     .domain([get2DMax(array2d), 0])
    //     // .domain([0,get2DMax(array2d)])
    //     // .interpolator(d3.scaleSequential(["red", "yellow", "green"]));
    //     // .interpolator(d3.interpolateSpectral)
    //     // .interpolator(d3.interpolateTurbo);

    lines.selectAll(".dot")
        .data(death => death)
        .enter()
        .append("circle") // Uses the enter().append() method
        .attr("class", "dot") // Assign a class for styling
        .attr("cx", function (_, i) { return xScale(i) })
        .attr("cy", function (d) { return yScale(d) })
        .attr("r", 2)
        // .style("stroke", (_, i) => lineColors(i))
        .style("stroke", "#333")
        .style("fill", "rgba(0,0,0,0.1)")

        .on("mouseover", function (event, d) {
            d3.select(this)
                .style("opacity", 0.7)

            tooltip.transition()
                .duration(200)
                .style("opacity", 1);
            tooltip.html(Math.round(d * 1000) / 1000)
                .style("left", (event.pageX) + "px")
                .style("top", (event.pageY - 28) + "px")
                .style("color", "white")
        })
        .on("mouseout", function (event, d) {
            d3.select(this)
                .style("opacity", 1)
            tooltip.transition()
                .duration(500)
                .style("opacity", 0);
        });


    // frame.append("path")
    //     // .datum(deaths[1])
    //     .datum(deaths[1])
    //     .attr("class", "line")
    //     .attr("d", line)
    // .line {
    //     fill: none;
    //     stroke: #ffab00;
    //     stroke-width: 3;
    // }

}

function barplot(frame, arrays, rwid) {

    let maxBarH = d3.max(arrays, array => d3.sum(array));
    console.log(maxBarH);
    // let data = Object.values(arrays)
    console.log(arrays);

    // arrays = arrays[age][cod]

    let data = d3.layout.stack()(arrays);
    console.log(data);

    let colorant = d3.scaleSequential()
        .domain([0, arrays.length])
        .interpolator(d3.interpolateSpectral)

    frame
        .selectAll("g")
        .data(arrays)
        .enter()
        .append("g")
        .attr("transform", (_, i) => `translate(${rwid * i},0)`)
        .selectAll("rect")
        .data(array => array)
        .enter()
        .append("rect")
        .attr("height", d => d/maxBarH * 1000)
        .attr("kay", _ => console.log(_))
        .attr("width", rwid - 0.3)
        .attr("fill", (_, i) => colorant(i))

}
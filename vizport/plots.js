const zip = (arr1, arr2) => arr1.map((k, i) => [k, arr2[i]]);


function heatmap(frame, data, rwid, scale) {

    frame
        .selectAll("g")
        .data(data)
        .join("g")
        .attr("transform", (_, i) => `translate(0, ${i * rwid})`)

        .selectAll("rect")
        .data(d => d)
        .join("rect")
        .attr("width", rwid - 0.3) // -0.5
        .attr("height", rwid - 0.3) // -0.5
        .attr("x", (_, i) => rwid * i)
        .style("fill", d => scale(d))
        // .attr("order", (_, i) => i)
        // .on("mouseover", function () {
        //     let order = d3.select(this).attr("order");
        //     d3.select("span#locus").text(order);
        // })

}

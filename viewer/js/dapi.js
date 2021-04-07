function dapi(cfg) {
    console.log('Doing Dapi plot');

    var img = cfg.imageSize,
        map_tiles = cfg.map_tiles,
        roi = cfg.roi;

    var a = img[0] / (roi.x1 - roi.x0),
        b = -img[0] / (roi.x1 - roi.x0) * roi.x0,
        c = img[1] / (roi.y1 - roi.y0),
        d = -img[1] / (roi.y1 - roi.y0) * roi.y0;

    // This transformation maps a point from the roi domain to the domain defined by [0,0] amd [img[0], img[1]].
    var t = new L.Transformation(a, b, c, d);

    // The transformation in this CRS maps the the top left corner to (0,0) and the bottom right to (256, 256)
    L.CRS.MySimple = L.extend({}, L.CRS.Simple, {
        transformation: new L.Transformation(1 / 1024, 0, 1 / 1024, 0),
    });


    var minZoom = 0,
        maxZoom = 8;
    var dapiLayer = L.tileLayer(map_tiles, {minZoom: minZoom, maxZoom: maxZoom});

    map = L.map('mymap', {
        layers: [dapiLayer],
        crs: L.CRS.MySimple,
        attributionControl: false,
    }).setView([img[1], img[0] / 2], 2);


    var dapiData = {};
    dapiData.map = map;
    dapiData.t = t;
    dapiData.minZoom = minZoom;
    dapiData.maxZoom = maxZoom;

    return dapiData
}


function dapiChart(config) {

    dapiConfig = dapi(config);
    var map = dapiConfig.map,
        minZoom = dapiConfig.minZoom,
        maxZoom = dapiConfig.maxZoom;

    // var gene_tiles = [];
    var gene_names;


    const repo = (gene) => {
        // return "https://raw.githubusercontent.com/acycliq/gene_tiles/master/" + gene + "/{z}/{y}/{x}.png"
        return "./src/pyramid/" + gene + "/{z}/{y}/{x}.png"
    };

    gene_names = glyphSettings().map(d => d.gene);
    // gene_names.forEach(gene => {
    //     var d = {};
    //     d.gene = gene;
    //     d.tileLayer = L.tileLayer(repo(gene),  {minZoom: minZoom, maxZoom: maxZoom});
    //     gene_tiles.push(d)
    // });

    var geneLayers_dict = {};
    gene_names.forEach(gene => {
        geneLayers_dict[gene] = L.tileLayer(repo(gene),  {minZoom: minZoom, maxZoom: maxZoom});
    });

    var Tph2Layer = L.tileLayer('./src/pyramid/Tph2/{z}/{y}/{x}.png', {minZoom: minZoom, maxZoom: maxZoom});
    var Slc1a2Layer = L.tileLayer('./src/pyramid/Slc1a2/{z}/{y}/{x}.png', {minZoom: minZoom, maxZoom: maxZoom});
    var Nos1Layer = L.tileLayer('./src/pyramid/Nos1/{z}/{y}/{x}.png', {minZoom: minZoom, maxZoom: maxZoom});

    // map.addLayer(Tph2Layer);
    // map.addLayer(Slc1a2Layer);
    // map.addLayer(Nos1Layer);

    gene_names.forEach(gene => {
        // map.addLayer(geneLayers_dict[gene])
    })

    // var baseLayers = {
    //     "Dapi image": dapiLayer,
    // };

    // var overlayMaps = {
    //     "Tph2": Tph2Layer,
    //     "Slc1a2": Slc1a2Layer,
    //     "Nos1": Nos1Layer,
    //
    // };

    // var overlayMaps = Object.assign({}, ...gene_tiles.map((x) => ({[x.gene]: x.tileLayer})));
    L.control.layers({}, geneLayers_dict, {collapsed: false}).addTo(map);

}


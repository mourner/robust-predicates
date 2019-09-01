
import {
    epsilon, splitter, resulterrbound, estimate,
    fast_expansion_sum_zeroelim, scale_expansion_zeroelim
} from './util.js';

const isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
const isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
const isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;

const ab = new Float64Array(4);
const bc = new Float64Array(4);
const cd = new Float64Array(4);
const de = new Float64Array(4);
const ea = new Float64Array(4);
const ac = new Float64Array(4);
const bd = new Float64Array(4);
const ce = new Float64Array(4);
const da = new Float64Array(4);
const eb = new Float64Array(4);
const temp8a = new Float64Array(8);
const temp8b = new Float64Array(8);
const temp16 = new Float64Array(16);
const abc = new Float64Array(24);
const bcd = new Float64Array(24);
const cde = new Float64Array(24);
const dea = new Float64Array(24);
const eab = new Float64Array(24);
const abd = new Float64Array(24);
const bce = new Float64Array(24);
const cda = new Float64Array(24);
const deb = new Float64Array(24);
const eac = new Float64Array(24);
const temp48a = new Float64Array(48);
const temp48b = new Float64Array(48);
const abcd = new Float64Array(96);
const bcde = new Float64Array(96);
const cdea = new Float64Array(96);
const deab = new Float64Array(96);
const eabc = new Float64Array(96);
const temp192 = new Float64Array(192);
const det384x = new Float64Array(384);
const det384y = new Float64Array(384);
const det384z = new Float64Array(384);
const detxy = new Float64Array(768);
const adet = new Float64Array(1152);
const bdet = new Float64Array(1152);
const cdet = new Float64Array(1152);
const ddet = new Float64Array(1152);
const edet = new Float64Array(1152);
const abdet = new Float64Array(2304);
const cddet = new Float64Array(2304);
const cdedet = new Float64Array(3456);
const deter = new Float64Array(5760);

function insphereexact(pa, pb, pc, pd, pe) {
    let axby1, bxcy1, cxdy1, dxey1, exay1;
    let bxay1, cxby1, dxcy1, exdy1, axey1;
    let axcy1, bxdy1, cxey1, dxay1, exby1;
    let cxay1, dxby1, excy1, axdy1, bxey1;
    let axby0, bxcy0, cxdy0, dxey0, exay0;
    let bxay0, cxby0, dxcy0, exdy0, axey0;
    let axcy0, bxdy0, cxey0, dxay0, exby0;
    let cxay0, dxby0, excy0, axdy0, bxey0;
    let temp8alen, temp8blen, temp16len;
    let temp48alen, temp48blen;
    let xlen, ylen, zlen;
    let xylen;
    let i;

    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0;

    $Two_Product(pa[0], pb[1], axby1, axby0);
    $Two_Product(pb[0], pa[1], bxay1, bxay0);
    $Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

    $Two_Product(pb[0], pc[1], bxcy1, bxcy0);
    $Two_Product(pc[0], pb[1], cxby1, cxby0);
    $Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

    $Two_Product(pc[0], pd[1], cxdy1, cxdy0);
    $Two_Product(pd[0], pc[1], dxcy1, dxcy0);
    $Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

    $Two_Product(pd[0], pe[1], dxey1, dxey0);
    $Two_Product(pe[0], pd[1], exdy1, exdy0);
    $Two_Two_Diff(dxey1, dxey0, exdy1, exdy0, de[3], de[2], de[1], de[0]);

    $Two_Product(pe[0], pa[1], exay1, exay0);
    $Two_Product(pa[0], pe[1], axey1, axey0);
    $Two_Two_Diff(exay1, exay0, axey1, axey0, ea[3], ea[2], ea[1], ea[0]);

    $Two_Product(pa[0], pc[1], axcy1, axcy0);
    $Two_Product(pc[0], pa[1], cxay1, cxay0);
    $Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

    $Two_Product(pb[0], pd[1], bxdy1, bxdy0);
    $Two_Product(pd[0], pb[1], dxby1, dxby0);
    $Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

    $Two_Product(pc[0], pe[1], cxey1, cxey0);
    $Two_Product(pe[0], pc[1], excy1, excy0);
    $Two_Two_Diff(cxey1, cxey0, excy1, excy0, ce[3], ce[2], ce[1], ce[0]);

    $Two_Product(pd[0], pa[1], dxay1, dxay0);
    $Two_Product(pa[0], pd[1], axdy1, axdy0);
    $Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

    $Two_Product(pe[0], pb[1], exby1, exby0);
    $Two_Product(pb[0], pe[1], bxey1, bxey0);
    $Two_Two_Diff(exby1, exby0, bxey1, bxey0, eb[3], eb[2], eb[1], eb[0]);

    temp8alen = scale_expansion_zeroelim(4, bc, pa[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, ac, -pb[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ab, pc[2], temp8a);
    const abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, abc);

    temp8alen = scale_expansion_zeroelim(4, cd, pb[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, bd, -pc[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, bc, pd[2], temp8a);
    const bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, bcd);

    temp8alen = scale_expansion_zeroelim(4, de, pc[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, ce, -pd[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, cd, pe[2], temp8a);
    let cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cde);

    temp8alen = scale_expansion_zeroelim(4, ea, pd[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, da, -pe[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, de, pa[2], temp8a);
    const dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, dea);

    temp8alen = scale_expansion_zeroelim(4, ab, pe[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, eb, -pa[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ea, pb[2], temp8a);
    const eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, eab);

    temp8alen = scale_expansion_zeroelim(4, bd, pa[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, da, pb[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ab, pd[2], temp8a);
    const abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, abd);

    temp8alen = scale_expansion_zeroelim(4, ce, pb[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, eb, pc[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, bc, pe[2], temp8a);
    const bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, bce);

    temp8alen = scale_expansion_zeroelim(4, da, pc[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, ac, pd[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, cd, pa[2], temp8a);
    const cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cda);

    temp8alen = scale_expansion_zeroelim(4, eb, pd[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, bd, pe[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, de, pb[2], temp8a);
    const deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, deb);

    temp8alen = scale_expansion_zeroelim(4, ac, pe[2], temp8a);
    temp8blen = scale_expansion_zeroelim(4, ce, pa[2], temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ea, pc[2], temp8a);
    const eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, eac);

    temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, bcde);
    xlen = scale_expansion_zeroelim(bcdelen, bcde, pa[0], temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, pa[0], det384x);
    ylen = scale_expansion_zeroelim(bcdelen, bcde, pa[1], temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, pa[1], det384y);
    zlen = scale_expansion_zeroelim(bcdelen, bcde, pa[2], temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, pa[2], det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const alen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, adet);

    temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, cdea);
    xlen = scale_expansion_zeroelim(cdealen, cdea, pb[0], temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, pb[0], det384x);
    ylen = scale_expansion_zeroelim(cdealen, cdea, pb[1], temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, pb[1], det384y);
    zlen = scale_expansion_zeroelim(cdealen, cdea, pb[2], temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, pb[2], det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const blen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, bdet);

    temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, deab);
    xlen = scale_expansion_zeroelim(deablen, deab, pc[0], temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, pc[0], det384x);
    ylen = scale_expansion_zeroelim(deablen, deab, pc[1], temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, pc[1], det384y);
    zlen = scale_expansion_zeroelim(deablen, deab, pc[2], temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, pc[2], det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const clen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cdet);

    temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, eabc);
    xlen = scale_expansion_zeroelim(eabclen, eabc, pd[0], temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, pd[0], det384x);
    ylen = scale_expansion_zeroelim(eabclen, eabc, pd[1], temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, pd[1], det384y);
    zlen = scale_expansion_zeroelim(eabclen, eabc, pd[2], temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, pd[2], det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const dlen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, ddet);

    temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, abcd);
    xlen = scale_expansion_zeroelim(abcdlen, abcd, pe[0], temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, pe[0], det384x);
    ylen = scale_expansion_zeroelim(abcdlen, abcd, pe[1], temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, pe[1], det384y);
    zlen = scale_expansion_zeroelim(abcdlen, abcd, pe[2], temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, pe[2], det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const elen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, edet);

    const ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    const cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
    cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cdedet);
    const deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, deter);

    return deter[deterlen - 1];
}

const temp8c = new Float64Array(8);
const temp24 = new Float64Array(24);
const temp48 = new Float64Array(48);
const xdet = new Float64Array(96);
const ydet = new Float64Array(96);
const zdet = new Float64Array(96);
const xydet = new Float64Array(192);
const fin1 = new Float64Array(1152);

function insphereadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez, permanent) {
    let aexbey1, bexaey1, bexcey1, cexbey1;
    let cexdey1, dexcey1, dexaey1, aexdey1;
    let aexcey1, cexaey1, bexdey1, dexbey1;
    let aexbey0, bexaey0, bexcey0, cexbey0;
    let cexdey0, dexcey0, dexaey0, aexdey0;
    let aexcey0, cexaey0, bexdey0, dexbey0;
    let ab3, bc3, cd3, da3, ac3, bd3;
    let temp8alen, temp8blen, temp8clen, temp16len, temp24len, temp48len;
    let xlen, ylen, zlen, xylen;

    let aextail, bextail, cextail, dextail;
    let aeytail, beytail, ceytail, deytail;
    let aeztail, beztail, ceztail, deztail;

    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0;

    const aex = ax - ex;
    const bex = bx - ex;
    const cex = cx - ex;
    const dex = dx - ex;
    const aey = ay - ey;
    const bey = by - ey;
    const cey = cy - ey;
    const dey = dy - ey;
    const aez = az - ez;
    const bez = bz - ez;
    const cez = cz - ez;
    const dez = dz - ez;

    $Two_Product(aex, bey, aexbey1, aexbey0);
    $Two_Product(bex, aey, bexaey1, bexaey0);
    $Two_Two_Diff(aexbey1, aexbey0, bexaey1, bexaey0, ab3, ab[2], ab[1], ab[0]);
    ab[3] = ab3;

    $Two_Product(bex, cey, bexcey1, bexcey0);
    $Two_Product(cex, bey, cexbey1, cexbey0);
    $Two_Two_Diff(bexcey1, bexcey0, cexbey1, cexbey0, bc3, bc[2], bc[1], bc[0]);
    bc[3] = bc3;

    $Two_Product(cex, dey, cexdey1, cexdey0);
    $Two_Product(dex, cey, dexcey1, dexcey0);
    $Two_Two_Diff(cexdey1, cexdey0, dexcey1, dexcey0, cd3, cd[2], cd[1], cd[0]);
    cd[3] = cd3;

    $Two_Product(dex, aey, dexaey1, dexaey0);
    $Two_Product(aex, dey, aexdey1, aexdey0);
    $Two_Two_Diff(dexaey1, dexaey0, aexdey1, aexdey0, da3, da[2], da[1], da[0]);
    da[3] = da3;

    $Two_Product(aex, cey, aexcey1, aexcey0);
    $Two_Product(cex, aey, cexaey1, cexaey0);
    $Two_Two_Diff(aexcey1, aexcey0, cexaey1, cexaey0, ac3, ac[2], ac[1], ac[0]);
    ac[3] = ac3;

    $Two_Product(bex, dey, bexdey1, bexdey0);
    $Two_Product(dex, bey, dexbey1, dexbey0);
    $Two_Two_Diff(bexdey1, bexdey0, dexbey1, dexbey0, bd3, bd[2], bd[1], bd[0]);
    bd[3] = bd3;

    temp8alen = scale_expansion_zeroelim(4, cd, bez, temp8a);
    temp8blen = scale_expansion_zeroelim(4, bd, -cez, temp8b);
    temp8clen = scale_expansion_zeroelim(4, bc, dez, temp8c);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, temp24);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, aex, temp48);
    xlen = scale_expansion_zeroelim(temp48len, temp48, -aex, xdet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, aey, temp48);
    ylen = scale_expansion_zeroelim(temp48len, temp48, -aey, ydet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, aez, temp48);
    zlen = scale_expansion_zeroelim(temp48len, temp48, -aez, zdet);
    xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
    const alen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, adet);

    temp8alen = scale_expansion_zeroelim(4, da, cez, temp8a);
    temp8blen = scale_expansion_zeroelim(4, ac, dez, temp8b);
    temp8clen = scale_expansion_zeroelim(4, cd, aez, temp8c);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, temp24);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, bex, temp48);
    xlen = scale_expansion_zeroelim(temp48len, temp48, bex, xdet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, bey, temp48);
    ylen = scale_expansion_zeroelim(temp48len, temp48, bey, ydet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, bez, temp48);
    zlen = scale_expansion_zeroelim(temp48len, temp48, bez, zdet);
    xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
    const blen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, bdet);

    temp8alen = scale_expansion_zeroelim(4, ab, dez, temp8a);
    temp8blen = scale_expansion_zeroelim(4, bd, aez, temp8b);
    temp8clen = scale_expansion_zeroelim(4, da, bez, temp8c);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, temp24);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, cex, temp48);
    xlen = scale_expansion_zeroelim(temp48len, temp48, -cex, xdet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, cey, temp48);
    ylen = scale_expansion_zeroelim(temp48len, temp48, -cey, ydet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, cez, temp48);
    zlen = scale_expansion_zeroelim(temp48len, temp48, -cez, zdet);
    xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
    const clen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cdet);

    temp8alen = scale_expansion_zeroelim(4, bc, aez, temp8a);
    temp8blen = scale_expansion_zeroelim(4, ac, -bez, temp8b);
    temp8clen = scale_expansion_zeroelim(4, ab, cez, temp8c);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c, temp16len, temp16, temp24);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, dex, temp48);
    xlen = scale_expansion_zeroelim(temp48len, temp48, dex, xdet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, dey, temp48);
    ylen = scale_expansion_zeroelim(temp48len, temp48, dey, ydet);
    temp48len = scale_expansion_zeroelim(temp24len, temp24, dez, temp48);
    zlen = scale_expansion_zeroelim(temp48len, temp48, dez, zdet);
    xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
    const dlen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, ddet);

    const ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    const cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
    const finlength = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, fin1);

    let det = estimate(finlength, fin1);
    let errbound = isperrboundB * permanent;
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    $Two_Diff_Tail(ax, ex, aex, aextail);
    $Two_Diff_Tail(ay, ey, aey, aeytail);
    $Two_Diff_Tail(az, ez, aez, aeztail);
    $Two_Diff_Tail(bx, ex, bex, bextail);
    $Two_Diff_Tail(by, ey, bey, beytail);
    $Two_Diff_Tail(bz, ez, bez, beztail);
    $Two_Diff_Tail(cx, ex, cex, cextail);
    $Two_Diff_Tail(cy, ey, cey, ceytail);
    $Two_Diff_Tail(cz, ez, cez, ceztail);
    $Two_Diff_Tail(dx, ex, dex, dextail);
    $Two_Diff_Tail(dy, ey, dey, deytail);
    $Two_Diff_Tail(dz, ez, dez, deztail);
    if ((aextail === 0.0) && (aeytail === 0.0) && (aeztail === 0.0) &&
        (bextail === 0.0) && (beytail === 0.0) && (beztail === 0.0) &&
        (cextail === 0.0) && (ceytail === 0.0) && (ceztail === 0.0) &&
        (dextail === 0.0) && (deytail === 0.0) && (deztail === 0.0)) {
        return det;
    }

    errbound = isperrboundC * permanent + resulterrbound * Math.abs(det);

    const abeps = (aex * beytail + bey * aextail) - (aey * bextail + bex * aeytail);
    const bceps = (bex * ceytail + cey * bextail) - (bey * cextail + cex * beytail);
    const cdeps = (cex * deytail + dey * cextail) - (cey * dextail + dex * ceytail);
    const daeps = (dex * aeytail + aey * dextail) - (dey * aextail + aex * deytail);
    const aceps = (aex * ceytail + cey * aextail) - (aey * cextail + cex * aeytail);
    const bdeps = (bex * deytail + dey * bextail) - (bey * dextail + dex * beytail);
    det +=
        (((bex * bex + bey * bey + bez * bez) * ((cez * daeps + dez * aceps + aez * cdeps) +
        (ceztail * da3 + deztail * ac3 + aeztail * cd3)) + (dex * dex + dey * dey + dez * dez) *
        ((aez * bceps - bez * aceps + cez * abeps) + (aeztail * bc3 - beztail * ac3 + ceztail * ab3))) -
        ((aex * aex + aey * aey + aez * aez) * ((bez * cdeps - cez * bdeps + dez * bceps) +
        (beztail * cd3 - ceztail * bd3 + deztail * bc3)) + (cex * cex + cey * cey + cez * cez) *
        ((dez * abeps + aez * bdeps + bez * daeps) + (deztail * ab3 + aeztail * bd3 + beztail * da3)))) +
        2.0 * (((bex * bextail + bey * beytail + bez * beztail) * (cez * da3 + dez * ac3 + aez * cd3) +
        (dex * dextail + dey * deytail + dez * deztail) * (aez * bc3 - bez * ac3 + cez * ab3)) -
        ((aex * aextail + aey * aeytail + aez * aeztail) * (bez * cd3 - cez * bd3 + dez * bc3) +
        (cex * cextail + cey * ceytail + cez * ceztail) * (dez * ab3 + aez * bd3 + bez * da3)));

    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    return insphereexact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez);
}

export function insphere(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez) {
    const aex = ax - ex;
    const bex = bx - ex;
    const cex = cx - ex;
    const dex = dx - ex;
    const aey = ay - ey;
    const bey = by - ey;
    const cey = cy - ey;
    const dey = dy - ey;
    const aez = az - ez;
    const bez = bz - ez;
    const cez = cz - ez;
    const dez = dz - ez;

    const aexbey = aex * bey;
    const bexaey = bex * aey;
    const ab = aexbey - bexaey;
    const bexcey = bex * cey;
    const cexbey = cex * bey;
    const bc = bexcey - cexbey;
    const cexdey = cex * dey;
    const dexcey = dex * cey;
    const cd = cexdey - dexcey;
    const dexaey = dex * aey;
    const aexdey = aex * dey;
    const da = dexaey - aexdey;

    const aexcey = aex * cey;
    const cexaey = cex * aey;
    const ac = aexcey - cexaey;
    const bexdey = bex * dey;
    const dexbey = dex * bey;
    const bd = bexdey - dexbey;

    const abc = aez * bc - bez * ac + cez * ab;
    const bcd = bez * cd - cez * bd + dez * bc;
    const cda = cez * da + dez * ac + aez * cd;
    const dab = dez * ab + aez * bd + bez * da;

    const alift = aex * aex + aey * aey + aez * aez;
    const blift = bex * bex + bey * bey + bez * bez;
    const clift = cex * cex + cey * cey + cez * cez;
    const dlift = dex * dex + dey * dey + dez * dez;

    const det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    const aezplus = Math.abs(aez);
    const bezplus = Math.abs(bez);
    const cezplus = Math.abs(cez);
    const dezplus = Math.abs(dez);
    const aexbeyplus = Math.abs(aexbey);
    const bexaeyplus = Math.abs(bexaey);
    const bexceyplus = Math.abs(bexcey);
    const cexbeyplus = Math.abs(cexbey);
    const cexdeyplus = Math.abs(cexdey);
    const dexceyplus = Math.abs(dexcey);
    const dexaeyplus = Math.abs(dexaey);
    const aexdeyplus = Math.abs(aexdey);
    const aexceyplus = Math.abs(aexcey);
    const cexaeyplus = Math.abs(cexaey);
    const bexdeyplus = Math.abs(bexdey);
    const dexbeyplus = Math.abs(dexbey);
    const permanent =
        ((cexdeyplus + dexceyplus) * bezplus + (dexbeyplus + bexdeyplus) * cezplus + (bexceyplus + cexbeyplus) * dezplus) * alift +
        ((dexaeyplus + aexdeyplus) * cezplus + (aexceyplus + cexaeyplus) * dezplus + (cexdeyplus + dexceyplus) * aezplus) * blift +
        ((aexbeyplus + bexaeyplus) * dezplus + (bexdeyplus + dexbeyplus) * aezplus + (dexaeyplus + aexdeyplus) * bezplus) * clift +
        ((bexceyplus + cexbeyplus) * aezplus + (cexaeyplus + aexceyplus) * bezplus + (aexbeyplus + bexaeyplus) * cezplus) * dlift;
    const errbound = isperrboundA * permanent;
    if ((det > errbound) || (-det > errbound)) {
        return det;
    }

    return insphereadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez, permanent);
}

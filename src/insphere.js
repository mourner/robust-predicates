import {
    epsilon, splitter, resulterrbound, estimate, vec,
    fast_expansion_sum_zeroelim, scale_expansion_zeroelim
} from './util.js';

const isperrboundA = (16 + 224 * epsilon) * epsilon;
const isperrboundB = (5 + 72 * epsilon) * epsilon;
const isperrboundC = (71 + 1408 * epsilon) * epsilon * epsilon;

const ab = vec(4);
const bc = vec(4);
const cd = vec(4);
const de = vec(4);
const ea = vec(4);
const ac = vec(4);
const bd = vec(4);
const ce = vec(4);
const da = vec(4);
const eb = vec(4);
const temp8a = vec(8);
const temp8b = vec(8);
const temp16 = vec(16);
const abc = vec(24);
const bcd = vec(24);
const cde = vec(24);
const dea = vec(24);
const eab = vec(24);
const abd = vec(24);
const bce = vec(24);
const cda = vec(24);
const deb = vec(24);
const eac = vec(24);
const temp48a = vec(48);
const temp48b = vec(48);
const abcd = vec(96);
const bcde = vec(96);
const cdea = vec(96);
const deab = vec(96);
const eabc = vec(96);
const temp192 = vec(192);
const det384x = vec(384);
const det384y = vec(384);
const det384z = vec(384);
const detxy = vec(768);
const adet = vec(1152);
const bdet = vec(1152);
const cdet = vec(1152);
const ddet = vec(1152);
const edet = vec(1152);
const abdet = vec(2304);
const cddet = vec(2304);
const cdedet = vec(3456);
const deter = vec(5760);

function insphereexact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez) {
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;
    let temp8alen, temp8blen, temp16len;
    let temp48alen, temp48blen;
    let xlen, ylen, zlen, xylen, i;

    $Cross_Product(ax, ay, bx, by, ab);
    $Cross_Product(bx, by, cx, cy, bc);
    $Cross_Product(cx, cy, dx, dy, cd);
    $Cross_Product(dx, dy, ex, ey, de);
    $Cross_Product(ex, ey, ax, ay, ea);
    $Cross_Product(ax, ay, cx, cy, ac);
    $Cross_Product(bx, by, dx, dy, bd);
    $Cross_Product(cx, cy, ex, ey, ce);
    $Cross_Product(dx, dy, ax, ay, da);
    $Cross_Product(ex, ey, bx, by, eb);

    temp8alen = scale_expansion_zeroelim(4, bc, az, temp8a);
    temp8blen = scale_expansion_zeroelim(4, ac, -bz, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ab, cz, temp8a);
    const abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, abc);

    temp8alen = scale_expansion_zeroelim(4, cd, bz, temp8a);
    temp8blen = scale_expansion_zeroelim(4, bd, -cz, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, bc, dz, temp8a);
    const bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, bcd);

    temp8alen = scale_expansion_zeroelim(4, de, cz, temp8a);
    temp8blen = scale_expansion_zeroelim(4, ce, -dz, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, cd, ez, temp8a);
    let cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cde);

    temp8alen = scale_expansion_zeroelim(4, ea, dz, temp8a);
    temp8blen = scale_expansion_zeroelim(4, da, -ez, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, de, az, temp8a);
    const dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, dea);

    temp8alen = scale_expansion_zeroelim(4, ab, ez, temp8a);
    temp8blen = scale_expansion_zeroelim(4, eb, -az, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ea, bz, temp8a);
    const eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, eab);

    temp8alen = scale_expansion_zeroelim(4, bd, az, temp8a);
    temp8blen = scale_expansion_zeroelim(4, da, bz, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ab, dz, temp8a);
    const abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, abd);

    temp8alen = scale_expansion_zeroelim(4, ce, bz, temp8a);
    temp8blen = scale_expansion_zeroelim(4, eb, cz, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, bc, ez, temp8a);
    const bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, bce);

    temp8alen = scale_expansion_zeroelim(4, da, cz, temp8a);
    temp8blen = scale_expansion_zeroelim(4, ac, dz, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, cd, az, temp8a);
    const cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, cda);

    temp8alen = scale_expansion_zeroelim(4, eb, dz, temp8a);
    temp8blen = scale_expansion_zeroelim(4, bd, ez, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, de, bz, temp8a);
    const deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, deb);

    temp8alen = scale_expansion_zeroelim(4, ac, ez, temp8a);
    temp8blen = scale_expansion_zeroelim(4, ce, az, temp8b);
    temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b, temp16);
    temp8alen = scale_expansion_zeroelim(4, ea, cz, temp8a);
    const eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16, eac);

    temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, bcde);
    xlen = scale_expansion_zeroelim(bcdelen, bcde, ax, temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, ax, det384x);
    ylen = scale_expansion_zeroelim(bcdelen, bcde, ay, temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, ay, det384y);
    zlen = scale_expansion_zeroelim(bcdelen, bcde, az, temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, az, det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const alen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, adet);

    temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, cdea);
    xlen = scale_expansion_zeroelim(cdealen, cdea, bx, temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, bx, det384x);
    ylen = scale_expansion_zeroelim(cdealen, cdea, by, temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, by, det384y);
    zlen = scale_expansion_zeroelim(cdealen, cdea, bz, temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, bz, det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const blen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, bdet);

    temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, deab);
    xlen = scale_expansion_zeroelim(deablen, deab, cx, temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, cx, det384x);
    ylen = scale_expansion_zeroelim(deablen, deab, cy, temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, cy, det384y);
    zlen = scale_expansion_zeroelim(deablen, deab, cz, temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, cz, det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const clen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cdet);

    temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, eabc);
    xlen = scale_expansion_zeroelim(eabclen, eabc, dx, temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, dx, det384x);
    ylen = scale_expansion_zeroelim(eabclen, eabc, dy, temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, dy, det384y);
    zlen = scale_expansion_zeroelim(eabclen, eabc, dz, temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, dz, det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const dlen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, ddet);

    temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a);
    temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b);
    for (i = 0; i < temp48blen; i++) {
        temp48b[i] = -temp48b[i];
    }
    const abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a, temp48blen, temp48b, abcd);
    xlen = scale_expansion_zeroelim(abcdlen, abcd, ex, temp192);
    xlen = scale_expansion_zeroelim(xlen, temp192, ex, det384x);
    ylen = scale_expansion_zeroelim(abcdlen, abcd, ey, temp192);
    ylen = scale_expansion_zeroelim(ylen, temp192, ey, det384y);
    zlen = scale_expansion_zeroelim(abcdlen, abcd, ez, temp192);
    zlen = scale_expansion_zeroelim(zlen, temp192, ez, det384z);
    xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
    const elen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, edet);

    const ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    const cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
    cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cdedet);
    const deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, deter);

    return deter[deterlen - 1];
}

const temp8c = vec(8);
const temp24 = vec(24);
const temp48 = vec(48);
const xdet = vec(96);
const ydet = vec(96);
const zdet = vec(96);
const xydet = vec(192);
const fin1 = vec(1152);

function insphereadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez, permanent) {
    let ab3, bc3, cd3, da3, ac3, bd3;
    let temp8alen, temp8blen, temp8clen, temp16len, temp24len, temp48len;
    let xlen, ylen, zlen, xylen;

    let aextail, bextail, cextail, dextail;
    let aeytail, beytail, ceytail, deytail;
    let aeztail, beztail, ceztail, deztail;

    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0;

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

    $Cross_Product(aex, aey, bex, bey, ab, ab3);
    $Cross_Product(bex, bey, cex, cey, bc, bc3);
    $Cross_Product(cex, cey, dex, dey, cd, cd3);
    $Cross_Product(dex, dey, aex, aey, da, da3);
    $Cross_Product(aex, aey, cex, cey, ac, ac3);
    $Cross_Product(bex, bey, dex, dey, bd, bd3);

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
    if (det >= errbound || -det >= errbound) {
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
    if (aextail === 0 && aeytail === 0 && aeztail === 0 &&
        bextail === 0 && beytail === 0 && beztail === 0 &&
        cextail === 0 && ceytail === 0 && ceztail === 0 &&
        dextail === 0 && deytail === 0 && deztail === 0) {
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
        2 * (((bex * bextail + bey * beytail + bez * beztail) * (cez * da3 + dez * ac3 + aez * cd3) +
        (dex * dextail + dey * deytail + dez * deztail) * (aez * bc3 - bez * ac3 + cez * ab3)) -
        ((aex * aextail + aey * aeytail + aez * aeztail) * (bez * cd3 - cez * bd3 + dez * bc3) +
        (cex * cextail + cey * ceytail + cez * ceztail) * (dez * ab3 + aez * bd3 + bez * da3)));

    if (det >= errbound || -det >= errbound) {
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

    const det = (clift * dab - dlift * abc) + (alift * bcd - blift * cda);

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
    if (det > errbound || -det > errbound) {
        return det;
    }
    return -insphereadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez, permanent);
}

export function inspherefast(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, pex, pey, pez) {
    const aex = pax - pex;
    const bex = pbx - pex;
    const cex = pcx - pex;
    const dex = pdx - pex;
    const aey = pay - pey;
    const bey = pby - pey;
    const cey = pcy - pey;
    const dey = pdy - pey;
    const aez = paz - pez;
    const bez = pbz - pez;
    const cez = pcz - pez;
    const dez = pdz - pez;

    const ab = aex * bey - bex * aey;
    const bc = bex * cey - cex * bey;
    const cd = cex * dey - dex * cey;
    const da = dex * aey - aex * dey;

    const ac = aex * cey - cex * aey;
    const bd = bex * dey - dex * bey;

    const abc = aez * bc - bez * ac + cez * ab;
    const bcd = bez * cd - cez * bd + dez * bc;
    const cda = cez * da + dez * ac + aez * cd;
    const dab = dez * ab + aez * bd + bez * da;

    const alift = aex * aex + aey * aey + aez * aez;
    const blift = bex * bex + bey * bey + bez * bez;
    const clift = cex * cex + cey * cey + cez * cez;
    const dlift = dex * dex + dey * dey + dez * dez;

    return (clift * dab - dlift * abc) + (alift * bcd - blift * cda);
}

import {
    epsilon, splitter, resulterrbound, estimate, vec,
    expansion_sum as sum, scale_expansion as scale
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
const abcd = vec(96);
const bcde = vec(96);
const cdea = vec(96);
const deab = vec(96);
const eabc = vec(96);
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

const _8 = vec(8);
const _8b = vec(8);
const _16 = vec(16);
const _24 = vec(24);
const _48 = vec(48);
const _48b = vec(48);
const _192 = vec(192);
const _384x = vec(384);
const _384y = vec(384);
const _384z = vec(384);

function insphereexact(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez) {
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;
    let _8len, _8blen, _16len, _48len, _48blen;
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

    _8len = scale(4, bc, az, _8);
    _8blen = scale(4, ac, -bz, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, ab, cz, _8);
    const abclen = sum(_8len, _8, _16len, _16, abc);

    _8len = scale(4, cd, bz, _8);
    _8blen = scale(4, bd, -cz, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, bc, dz, _8);
    const bcdlen = sum(_8len, _8, _16len, _16, bcd);

    _8len = scale(4, de, cz, _8);
    _8blen = scale(4, ce, -dz, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, cd, ez, _8);
    let cdelen = sum(_8len, _8, _16len, _16, cde);

    _8len = scale(4, ea, dz, _8);
    _8blen = scale(4, da, -ez, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, de, az, _8);
    const dealen = sum(_8len, _8, _16len, _16, dea);

    _8len = scale(4, ab, ez, _8);
    _8blen = scale(4, eb, -az, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, ea, bz, _8);
    const eablen = sum(_8len, _8, _16len, _16, eab);

    _8len = scale(4, bd, az, _8);
    _8blen = scale(4, da, bz, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, ab, dz, _8);
    const abdlen = sum(_8len, _8, _16len, _16, abd);

    _8len = scale(4, ce, bz, _8);
    _8blen = scale(4, eb, cz, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, bc, ez, _8);
    const bcelen = sum(_8len, _8, _16len, _16, bce);

    _8len = scale(4, da, cz, _8);
    _8blen = scale(4, ac, dz, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, cd, az, _8);
    const cdalen = sum(_8len, _8, _16len, _16, cda);

    _8len = scale(4, eb, dz, _8);
    _8blen = scale(4, bd, ez, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, de, bz, _8);
    const deblen = sum(_8len, _8, _16len, _16, deb);

    _8len = scale(4, ac, ez, _8);
    _8blen = scale(4, ce, az, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, ea, cz, _8);
    const eaclen = sum(_8len, _8, _16len, _16, eac);

    _48len = sum(cdelen, cde, bcelen, bce, _48);
    _48blen = sum(deblen, deb, bcdlen, bcd, _48b);
    for (i = 0; i < _48blen; i++) _48b[i] = -_48b[i];

    const bcdelen = sum(_48len, _48, _48blen, _48b, bcde);
    xlen = scale(scale(bcdelen, bcde, ax, _192), _192, ax, _384x);
    ylen = scale(scale(bcdelen, bcde, ay, _192), _192, ay, _384y);
    zlen = scale(scale(bcdelen, bcde, az, _192), _192, az, _384z);
    xylen = sum(xlen, _384x, ylen, _384y, detxy);
    const alen = sum(xylen, detxy, zlen, _384z, adet);

    _48len = sum(dealen, dea, cdalen, cda, _48);
    _48blen = sum(eaclen, eac, cdelen, cde, _48b);
    for (i = 0; i < _48blen; i++) _48b[i] = -_48b[i];

    const cdealen = sum(_48len, _48, _48blen, _48b, cdea);
    xlen = scale(scale(cdealen, cdea, bx, _192), _192, bx, _384x);
    ylen = scale(scale(cdealen, cdea, by, _192), _192, by, _384y);
    zlen = scale(scale(cdealen, cdea, bz, _192), _192, bz, _384z);
    xylen = sum(xlen, _384x, ylen, _384y, detxy);
    const blen = sum(xylen, detxy, zlen, _384z, bdet);

    _48len = sum(eablen, eab, deblen, deb, _48);
    _48blen = sum(abdlen, abd, dealen, dea, _48b);
    for (i = 0; i < _48blen; i++) _48b[i] = -_48b[i];

    const deablen = sum(_48len, _48, _48blen, _48b, deab);
    xlen = scale(scale(deablen, deab, cx, _192), _192, cx, _384x);
    ylen = scale(scale(deablen, deab, cy, _192), _192, cy, _384y);
    zlen = scale(scale(deablen, deab, cz, _192), _192, cz, _384z);
    xylen = sum(xlen, _384x, ylen, _384y, detxy);
    const clen = sum(xylen, detxy, zlen, _384z, cdet);

    _48len = sum(abclen, abc, eaclen, eac, _48);
    _48blen = sum(bcelen, bce, eablen, eab, _48b);
    for (i = 0; i < _48blen; i++) _48b[i] = -_48b[i];

    const eabclen = sum(_48len, _48, _48blen, _48b, eabc);
    xlen = scale(scale(eabclen, eabc, dx, _192), _192, dx, _384x);
    ylen = scale(scale(eabclen, eabc, dy, _192), _192, dy, _384y);
    zlen = scale(scale(eabclen, eabc, dz, _192), _192, dz, _384z);
    xylen = sum(xlen, _384x, ylen, _384y, detxy);
    const dlen = sum(xylen, detxy, zlen, _384z, ddet);

    _48len = sum(bcdlen, bcd, abdlen, abd, _48);
    _48blen = sum(cdalen, cda, abclen, abc, _48b);
    for (i = 0; i < _48blen; i++) _48b[i] = -_48b[i];

    const abcdlen = sum(_48len, _48, _48blen, _48b, abcd);
    xlen = scale(scale(abcdlen, abcd, ex, _192), _192, ex, _384x);
    ylen = scale(scale(abcdlen, abcd, ey, _192), _192, ey, _384y);
    zlen = scale(scale(abcdlen, abcd, ez, _192), _192, ez, _384z);
    xylen = sum(xlen, _384x, ylen, _384y, detxy);
    const elen = sum(xylen, detxy, zlen, _384z, edet);

    const ablen = sum(alen, adet, blen, bdet, abdet);
    const cdlen = sum(clen, cdet, dlen, ddet, cddet);
    cdelen = sum(cdlen, cddet, elen, edet, cdedet);
    const deterlen = sum(ablen, abdet, cdelen, cdedet, deter);

    return deter[deterlen - 1];
}

const xdet = vec(96);
const ydet = vec(96);
const zdet = vec(96);
const xydet = vec(192);
const fin = vec(1152);

function insphereadapt(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, ex, ey, ez, permanent) {
    let ab3, bc3, cd3, da3, ac3, bd3;
    let _8len, _8blen, _16len, _24len, _48len;
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

    _8len = scale(4, cd, bez, _8);
    _8blen = scale(4, bd, -cez, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, bc, dez, _8);
    _24len = sum(_8len, _8, _16len, _16, _24);
    _48len = scale(_24len, _24, aex, _48);
    xlen = scale(_48len, _48, -aex, xdet);
    _48len = scale(_24len, _24, aey, _48);
    ylen = scale(_48len, _48, -aey, ydet);
    _48len = scale(_24len, _24, aez, _48);
    zlen = scale(_48len, _48, -aez, zdet);
    xylen = sum(xlen, xdet, ylen, ydet, xydet);
    const alen = sum(xylen, xydet, zlen, zdet, adet);

    _8len = scale(4, da, cez, _8);
    _8blen = scale(4, ac, dez, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, cd, aez, _8);
    _24len = sum(_8len, _8, _16len, _16, _24);
    _48len = scale(_24len, _24, bex, _48);
    xlen = scale(_48len, _48, bex, xdet);
    _48len = scale(_24len, _24, bey, _48);
    ylen = scale(_48len, _48, bey, ydet);
    _48len = scale(_24len, _24, bez, _48);
    zlen = scale(_48len, _48, bez, zdet);
    xylen = sum(xlen, xdet, ylen, ydet, xydet);
    const blen = sum(xylen, xydet, zlen, zdet, bdet);

    _8len = scale(4, ab, dez, _8);
    _8blen = scale(4, bd, aez, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, da, bez, _8);
    _24len = sum(_8len, _8, _16len, _16, _24);
    _48len = scale(_24len, _24, cex, _48);
    xlen = scale(_48len, _48, -cex, xdet);
    _48len = scale(_24len, _24, cey, _48);
    ylen = scale(_48len, _48, -cey, ydet);
    _48len = scale(_24len, _24, cez, _48);
    zlen = scale(_48len, _48, -cez, zdet);
    xylen = sum(xlen, xdet, ylen, ydet, xydet);
    const clen = sum(xylen, xydet, zlen, zdet, cdet);

    _8len = scale(4, bc, aez, _8);
    _8blen = scale(4, ac, -bez, _8b);
    _16len = sum(_8len, _8, _8blen, _8b, _16);
    _8len = scale(4, ab, cez, _8);
    _24len = sum(_8len, _8, _16len, _16, _24);
    _48len = scale(_24len, _24, dex, _48);
    xlen = scale(_48len, _48, dex, xdet);
    _48len = scale(_24len, _24, dey, _48);
    ylen = scale(_48len, _48, dey, ydet);
    _48len = scale(_24len, _24, dez, _48);
    zlen = scale(_48len, _48, dez, zdet);
    xylen = sum(xlen, xdet, ylen, ydet, xydet);
    const dlen = sum(xylen, xydet, zlen, zdet, ddet);

    const ablen = sum(alen, adet, blen, bdet, abdet);
    const cdlen = sum(clen, cdet, dlen, ddet, cddet);
    const finlen = sum(ablen, abdet, cdlen, cddet, fin);

    let det = estimate(finlen, fin);
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

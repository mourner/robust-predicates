
const Fast_Two_Sum = (a, b, x, y) => `
    ${x} = ${a} + ${b};
    ${y} = ${b} - (${x} - ${a});
`.trim();

const Two_Sum = (a, b, x, y) => `
    ${x} = ${a} + ${b};
    bvirt = ${x} - ${a};
    ${y} = (${a} - (${x} - bvirt)) + (${b} - bvirt);
`.trim();

const Two_Diff_Tail = (a, b, x, y) => `
    bvirt = ${a} - ${x};
    ${y} = (${a} - (${x} + bvirt)) + (bvirt - ${b});
`.trim();

const Two_Diff = (a, b, x, y) => `
    ${x} = ${a} - ${b};
    ${Two_Diff_Tail(a, b, x, y)}
`.trim();

const Split = (a, ahi, alo) => `
    c = splitter * ${a};
    ${ahi} = c - (c - ${a});
    ${alo} = ${a} - ${ahi};
`.trim();

const Two_Product = (a, b, x, y) => `
    ${x} = ${a} * ${b};
    ${Split(a, 'ahi', 'alo')}
    ${Split(b, 'bhi', 'blo')}
    ${y} = alo * blo - (${x} - ahi * bhi - alo * bhi - ahi * blo);
`.trim();

const Two_Product_Presplit = (a, b, bhi, blo, x, y) => `
    ${x} = ${a} * ${b};
    ${Split(a, 'ahi', 'alo')}
    ${y} = alo * ${blo} - (${x} - ahi * ${bhi} - alo * ${bhi} - ahi * ${blo});
`.trim();

const Square = (a, x, y) => `
    ${x} = ${a} * ${a};
    ${Split(a, 'ahi', 'alo')}
    ${y} = alo * alo - (${x} - ahi * ahi - (ahi + ahi) * alo);
`.trim();

const Two_One_Sum = (a1, a0, b, x2, x1, x0) => `
    ${Two_Sum(a0, b, '_i', x0)}
    ${Two_Sum(a1, '_i', x2, x1)}
`.trim();

const Two_One_Diff = (a1, a0, b, x2, x1, x0) => `
    ${Two_Diff(a0, b, '_i', x0)}
    ${Two_Sum(a1, '_i', x2, x1)}
`.trim();

const Two_Two_Sum = (a1, a0, b1, b0, x3, x2, x1, x0) => `
    ${Two_One_Sum(a1, a0, b0, '_j', '_0', x0)}
    ${Two_One_Sum('_j', '_0', b1, x3, x2, x1)}
`.trim();

const Two_Two_Diff = (a1, a0, b1, b0, x3, x2, x1, x0) => `
    ${Two_One_Diff(a1, a0, b0, '_j', '_0', x0)}
    ${Two_One_Diff('_j', '_0', b1, x3, x2, x1)}
`.trim();

const Two_One_Product = (a1, a0, b, x3, x2, x1, x0) => `
    ${Split(b, 'bhi', 'blo')}
    ${Two_Product_Presplit(a0, b, 'bhi', 'blo', '_i', x0)}
    ${Two_Product_Presplit(a1, b, 'bhi', 'blo', '_j', '_0')}
    ${Two_Sum('_i', '_0', '_k', x1)}
    ${Fast_Two_Sum('_j', '_k', x3, x2)}
`.trim();

const src = `
const epsilon = 1.1102230246251565e-16;
const splitter = 134217729;

const resulterrbound = (3 + 8 * epsilon) * epsilon;
const ccwerrboundA = (3 + 16 * epsilon) * epsilon;
const ccwerrboundB = (2 + 12 * epsilon) * epsilon;
const ccwerrboundC = (9 + 64 * epsilon) * epsilon * epsilon;

function fast_expansion_sum_zeroelim(elen, e, flen, f, h) {
    let Q, Qnew, hh, bvirt;
    let enow = e[0];
    let fnow = f[0];
    let eindex = 0;
    let findex = 0;
    if ((fnow > enow) === (fnow > -enow)) {
        Q = enow;
        enow = e[++eindex];
    } else {
        Q = fnow;
        fnow = f[++findex];
    }
    let hindex = 0;
    if ((eindex < elen) && (findex < flen)) {
        if ((fnow > enow) === (fnow > -enow)) {
            ${Fast_Two_Sum('enow', 'Q', 'Qnew', 'hh')}
            enow = e[++eindex];
        } else {
            ${Fast_Two_Sum('fnow', 'Q', 'Qnew', 'hh')}
            fnow = f[++findex];
        }
        Q = Qnew;
        if (hh !== 0) {
            h[hindex++] = hh;
        }
        while ((eindex < elen) && (findex < flen)) {
            if ((fnow > enow) === (fnow > -enow)) {
                ${Two_Sum('Q', 'enow', 'Qnew', 'hh')}
                enow = e[++eindex];
            } else {
                ${Two_Sum('Q', 'fnow', 'Qnew', 'hh')}
                fnow = f[++findex];
            }
            Q = Qnew;
            if (hh !== 0) {
                h[hindex++] = hh;
            }
        }
    }
    while (eindex < elen) {
        ${Two_Sum('Q', 'enow', 'Qnew', 'hh')};
        enow = e[++eindex];
        Q = Qnew;
        if (hh !== 0) {
            h[hindex++] = hh;
        }
    }
    while (findex < flen) {
        ${Two_Sum('Q', 'fnow', 'Qnew', 'hh')};
        fnow = f[++findex];
        Q = Qnew;
        if (hh !== 0) {
            h[hindex++] = hh;
        }
    }
    if ((Q !== 0) || (hindex === 0)) {
        h[hindex++] = Q;
    }
    return hindex;
}

function scale_expansion_zeroelim(elen, e, b, h) {
  let Q, sum, hh, product1, product0;
  let bvirt, c, ahi, alo, bhi, blo;

  ${Split('b', 'bhi', 'blo')}
  ${Two_Product_Presplit('e[0]', 'b', 'bhi', 'blo', 'Q', 'hh')}
  let hindex = 0;
  if (hh !== 0) {
    h[hindex++] = hh;
  }
  for (let eindex = 1; eindex < elen; eindex++) {
    let enow = e[eindex];
    ${Two_Product_Presplit('enow', 'b', 'bhi', 'blo', 'product1', 'product0')}
    ${Two_Sum('Q', 'product0', 'sum', 'hh')}
    if (hh !== 0) {
      h[hindex++] = hh;
    }
    ${Fast_Two_Sum('product1', 'sum', 'Q', 'hh')}
    if (hh !== 0) {
      h[hindex++] = hh;
    }
  }
  if ((Q !== 0.0) || (hindex === 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

function estimate(elen, e) {
    let Q = e[0];
    for (let i = 1; i < elen; i++) Q += e[i];
    return Q;
}

const B = new Float64Array(4);
const C1 = new Float64Array(8);
const C2 = new Float64Array(12);
const D = new Float64Array(16);
const u = new Float64Array(4);

function orient2dadapt(ax, ay, bx, by, cx, cy, detsum) {
    let acxtail, acytail, bcxtail, bcytail;
    let detleft, detright, detlefttail, detrighttail;
    let B3, u3, s1, t1, s0, t0;
    let bvirt, c, ahi, alo, bhi, blo, _i, _j, _0;

    const acx = ax - cx;
    const bcx = bx - cx;
    const acy = ay - cy;
    const bcy = by - cy;

    ${Two_Product('acx', 'bcy', 'detleft', 'detlefttail')}
    ${Two_Product('acy', 'bcx', 'detright', 'detrighttail')}

    ${Two_Two_Diff('detleft', 'detlefttail', 'detright', 'detrighttail', 'B3', 'B[2]', 'B[1]', 'B[0]')}
    B[3] = B3;

    let det = estimate(4, B);
    let errbound = ccwerrboundB * detsum;
    if (det >= errbound || -det >= errbound) return det;

    ${Two_Diff_Tail('ax', 'cx', 'acx', 'acxtail')};
    ${Two_Diff_Tail('bx', 'cx', 'bcx', 'bcxtail')};
    ${Two_Diff_Tail('ay', 'cy', 'acy', 'acytail')};
    ${Two_Diff_Tail('by', 'cy', 'bcy', 'bcytail')};

    if ((acxtail === 0) && (acytail === 0) && (bcxtail === 0) && (bcytail === 0)) {
        return det;
    }

    errbound = ccwerrboundC * detsum + resulterrbound * Math.abs(det);
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if ((det >= errbound) || (-det >= errbound)) {
        return det;
    }

    ${Two_Product('acxtail', 'bcy', 's1', 's0')}
    ${Two_Product('acytail', 'bcx', 't1', 't0')}
    ${Two_Two_Diff('s1', 's0', 't1', 't0', 'u3', 'u[2]', 'u[1]', 'u[0]')}
    u[3] = u3;
    const C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

    ${Two_Product('acx', 'bcytail', 's1', 's0')}
    ${Two_Product('acy', 'bcxtail', 't1', 't0')}
    ${Two_Two_Diff('s1', 's0', 't1', 't0', 'u3', 'u[2]', 'u[1]', 'u[0]')}
    u[3] = u3;
    const C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

    ${Two_Product('acxtail', 'bcytail', 's1', 's0')}
    ${Two_Product('acytail', 'bcxtail', 't1', 't0')}
    ${Two_Two_Diff('s1', 's0', 't1', 't0', 'u3', 'u[2]', 'u[1]', 'u[0]')}
    u[3] = u3;
    const Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

    return(D[Dlength - 1]);
}

export function orient2d(ax, ay, bx, by, cx, cy) {
    const detleft = (ax - cx) * (by - cy);
    const detright = (ay - cy) * (bx - cx);
    const det = detleft - detright;
    let detsum;

    if (detleft > 0) {
        if (detright <= 0) return det;
        else detsum = detleft + detright;

    } else if (detleft < 0) {
        if (detright >= 0) return det;
        else detsum = -detleft - detright;

    } else return det;

    const errbound = ccwerrboundA * detsum;
    if (det >= errbound || -det >= errbound) return det;

    return orient2dadapt(ax, ay, bx, by, cx, cy, detsum);
}
`;

console.log(src);

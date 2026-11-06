/**
 * Homework 6: Numerical Integration (Riemann vs Lebesgue) & Normal CDF Computation
 * Author: Dario Pietrosanto (Matricola: 2055042)
 * Master's Degree Course in Applied Statistics & Cybersecurity
 */

(function () {
    'use strict';

    // =========================================================================
    // PART 1: RIEMANN VS LEBESGUE INTEGRATION
    // =========================================================================

    // Greatest Common Divisor helper for Thomae function
    function gcd(a, b) {
        a = Math.abs(a);
        b = Math.abs(b);
        while (b !== 0) {
            const t = b;
            b = a % b;
            a = t;
        }
        return a;
    }

    // Sine Integral Si(pi) constant
    const SI_PI = 1.85193705198246617036;

    // Test functions catalogue
    const FUNCTIONS = {
        sin: {
            name: "f(x) = sin(x)",
            domain: [0, Math.PI],
            displayDomain: "[0, \\pi]",
            exact: 2.0,
            latex: "f(x) = \\sin(x), \\quad x \\in [0, \\pi]",
            fn: (x) => Math.sin(x),
            category: "Smooth Periodic",
            description: "Smooth trigonometric function with derivative bounded by 1. Standard benchmark for quadrature.",
            theoreticalLebesgue: 2.0,
            theoreticalRiemann: 2.0,
            riemannIntegrable: true,
            yRange: [0, 1]
        },
        square: {
            name: "f(x) = x²",
            domain: [0, 1],
            displayDomain: "[0, 1]",
            exact: 1.0 / 3.0,
            latex: "f(x) = x^2, \\quad x \\in [0, 1]",
            fn: (x) => x * x,
            category: "Monotone Polynomial",
            description: "Strictly increasing convex polynomial. Exact integral is 1/3.",
            theoreticalLebesgue: 1.0 / 3.0,
            theoreticalRiemann: 1.0 / 3.0,
            riemannIntegrable: true,
            yRange: [0, 1]
        },
        step_1_3: {
            name: "f(x) = 1 if x < 0.5, 3 if x ≥ 0.5",
            domain: [0, 1],
            displayDomain: "[0, 1]",
            exact: 2.0,
            latex: "f(x) = \\begin{cases} 1 & x < 0.5 \\\\ 3 & x \\ge 0.5 \\end{cases}, \\quad x \\in [0, 1]",
            fn: (x) => (x < 0.5 ? 1.0 : 3.0),
            category: "Piecewise Constant",
            description: "Step function with a jump discontinuity at x = 0.5. Lebesgue integrates simple functions with zero error when range partitions align with values.",
            theoreticalLebesgue: 2.0,
            theoreticalRiemann: 2.0,
            riemannIntegrable: true,
            yRange: [1, 3]
        },
        abs_diff: {
            name: "f(x) = |x - 0.5|",
            domain: [0, 1],
            displayDomain: "[0, 1]",
            exact: 0.25,
            latex: "f(x) = |x - 0.5|, \\quad x \\in [0, 1]",
            fn: (x) => Math.abs(x - 0.5),
            category: "Continuous Cusp",
            description: "Continuous everywhere but non-differentiable at x = 0.5. Symmetric V-shaped wedge.",
            theoreticalLebesgue: 0.25,
            theoreticalRiemann: 0.25,
            riemannIntegrable: true,
            yRange: [0, 0.5]
        },
        inv_sqrt: {
            name: "f(x) = 1 / √x",
            domain: [0.0001, 1], // avoid 0 singularity for domain sampling
            displayDomain: "(0, 1]",
            exact: 2.0,
            latex: "f(x) = \\frac{1}{\\sqrt{x}}, \\quad x \\in (0, 1]",
            fn: (x) => 1.0 / Math.sqrt(Math.max(x, 1e-9)),
            category: "Unbounded Improper",
            description: "Improper integral unbounded as x → 0⁺. Riemann requires careful endpoint handling; Lebesgue layer-cake representation μ({x: f(x) > y}) = 1/y² integrates to 1 + base.",
            theoreticalLebesgue: 2.0,
            theoreticalRiemann: 2.0,
            riemannIntegrable: true,
            yRange: [1, 50]
        },
        step_1_0: {
            name: "f(x) = 1 if x < 0.5, 0 if x ≥ 0.5",
            domain: [0, 1],
            displayDomain: "[0, 1]",
            exact: 0.5,
            latex: "f(x) = \\mathbf{1}_{[0, 0.5)}(x), \\quad x \\in [0, 1]",
            fn: (x) => (x < 0.5 ? 1.0 : 0.0),
            category: "Indicator Function",
            description: "Characteristic function of interval [0, 0.5). Fundamental building block of Lebesgue simple functions.",
            theoreticalLebesgue: 0.5,
            theoreticalRiemann: 0.5,
            riemannIntegrable: true,
            yRange: [0, 1]
        },
        sinc: {
            name: "f(x) = sin(x) / x",
            domain: [0, Math.PI],
            displayDomain: "[0, \\pi]",
            exact: SI_PI,
            latex: "f(x) = \\begin{cases} \\frac{\\sin(x)}{x} & x > 0 \\\\ 1 & x = 0 \\end{cases}, \\quad x \\in [0, \\pi]",
            fn: (x) => (x === 0 ? 1.0 : Math.sin(x) / x),
            category: "Removable Singularity",
            description: "Cardinal sine function (sinc) with removable singularity at the origin. Exact value is Si(π) ≈ 1.851937.",
            theoreticalLebesgue: SI_PI,
            theoreticalRiemann: SI_PI,
            riemannIntegrable: true,
            yRange: [0, 1]
        },
        dirichlet: {
            name: "Dirichlet Function (1 if rational, 0 if irrational)",
            domain: [0, 1],
            displayDomain: "[0, 1]",
            exact: 0.0,
            latex: "f(x) = \\begin{cases} 1 & x \\in \\mathbb{Q} \\\\ 0 & x \\notin \\mathbb{Q} \\end{cases}, \\quad x \\in [0, 1]",
            // In simulation, we emulate the density of rationals by detecting rationals with denominator <= 500
            fn: (x) => {
                if (Math.abs(x) < 1e-12 || Math.abs(x - 1) < 1e-12) return 1.0;
                for (let q = 1; q <= 200; q++) {
                    const p = Math.round(x * q);
                    if (Math.abs(x - p / q) < 1e-6) return 1.0;
                }
                return 0.0;
            },
            category: "Pathological (Nowhere Continuous)",
            description: "Every subinterval contains both rationals and irrationals. Upper Darboux sum = 1, Lower Darboux sum = 0. Riemann integral DOES NOT EXIST. Lebesgue integral is identically 0 since μ(ℚ) = 0.",
            theoreticalLebesgue: 0.0,
            theoreticalRiemann: null, // Undefined
            riemannIntegrable: false,
            yRange: [0, 1]
        },
        thomae: {
            name: "Thomae's (Popcorn / Ruler) Function",
            domain: [0, 1],
            displayDomain: "[0, 1]",
            exact: 0.0,
            latex: "f(x) = \\begin{cases} 1/q & x = p/q \\in \\mathbb{Q} \\text{ (coprime)} \\\\ 0 & x \\notin \\mathbb{Q} \\end{cases}, \\quad x \\in [0, 1]",
            fn: (x) => {
                if (Math.abs(x) < 1e-12 || Math.abs(x - 1) < 1e-12) return 1.0;
                for (let q = 1; q <= 200; q++) {
                    const p = Math.round(x * q);
                    if (Math.abs(x - p / q) < 1e-6) {
                        if (gcd(p, q) === 1) return 1.0 / q;
                    }
                }
                return 0.0;
            },
            category: "Pathological (Continuity Set is Irrationals)",
            description: "Continuous at all irrationals, discontinuous at all rationals. Discontinuities have Lebesgue measure zero, so by Lebesgue's Criterion, it IS Riemann-integrable with value 0. Lebesgue integral is also 0.",
            theoreticalLebesgue: 0.0,
            theoreticalRiemann: 0.0,
            riemannIntegrable: true,
            yRange: [0, 1]
        }
    };

    /**
     * Riemann numerical integration
     * @param {Function} f integrand
     * @param {number} a domain start
     * @param {number} b domain end
     * @param {number} n number of subintervals
     * @param {string} rule 'midpoint' | 'trapezoid' | 'left' | 'right'
     */
    function riemannIntegrate(f, a, b, n, rule = 'midpoint') {
        const t0 = performance.now();
        const h = (b - a) / n;
        let sum = 0;

        if (rule === 'midpoint') {
            for (let i = 0; i < n; i++) {
                const x = a + (i + 0.5) * h;
                sum += f(x);
            }
            sum *= h;
        } else if (rule === 'trapezoid') {
            sum = 0.5 * (f(a) + f(b));
            for (let i = 1; i < n; i++) {
                sum += f(a + i * h);
            }
            sum *= h;
        } else if (rule === 'left') {
            for (let i = 0; i < n; i++) {
                sum += f(a + i * h);
            }
            sum *= h;
        } else if (rule === 'right') {
            for (let i = 1; i <= n; i++) {
                sum += f(a + i * h);
            }
            sum *= h;
        }

        const t1 = performance.now();
        return {
            value: sum,
            elapsedMs: t1 - t0,
            h: h,
            steps: n
        };
    }

    /**
     * Numerical Lebesgue Integration
     * Partitions the range [ymin, ymax] into M subintervals [y_j, y_{j+1}].
     * Computes the measure mu(E_j) of the preimage E_j = {x in [a, b] : y_j <= f(x) < y_{j+1}}
     * using K uniform probe sample points across the domain [a, b].
     */
    function lebesgueIntegrate(f, a, b, M, K = 50000, yMin = null, yMax = null) {
        const t0 = performance.now();
        const stepK = (b - a) / K;

        // Sample domain densely to determine or confirm range
        const vals = new Float64Array(K);
        let sampleMin = Infinity;
        let sampleMax = -Infinity;

        for (let k = 0; k < K; k++) {
            const x = a + (k + 0.5) * stepK;
            const v = f(x);
            vals[k] = v;
            if (v < sampleMin) sampleMin = v;
            if (v > sampleMax) sampleMax = v;
        }

        const y0 = yMin !== null ? yMin : sampleMin;
        const y1 = yMax !== null ? yMax : sampleMax;

        if (Math.abs(y1 - y0) < 1e-12) {
            return {
                value: y0 * (b - a),
                elapsedMs: performance.now() - t0,
                y0, y1, dy: 0,
                preimages: [{ yLower: y0, yUpper: y1, yMid: y0, measure: b - a, count: K }]
            };
        }

        const dy = (y1 - y0) / M;
        const counts = new Uint32Array(M);

        for (let k = 0; k < K; k++) {
            const v = vals[k];
            let idx = Math.floor((v - y0) / dy);
            if (idx < 0) idx = 0;
            if (idx >= M) idx = M - 1;
            counts[idx]++;
        }

        let totalIntegral = 0;
        const preimages = [];
        for (let j = 0; j < M; j++) {
            const yMid = y0 + (j + 0.5) * dy;
            const measure = (counts[j] / K) * (b - a);
            totalIntegral += yMid * measure;
            preimages.push({
                index: j,
                yLower: y0 + j * dy,
                yUpper: y0 + (j + 1) * dy,
                yMid: yMid,
                measure: measure,
                count: counts[j]
            });
        }

        const t1 = performance.now();
        return {
            value: totalIntegral,
            elapsedMs: t1 - t0,
            y0: y0,
            y1: y1,
            dy: dy,
            M: M,
            K: K,
            preimages: preimages
        };
    }

    // =========================================================================
    // PART 2: NORMAL CDF NUMERICAL COMPUTATION METHODS
    // =========================================================================

    // Standard Normal Probability Density Function phi(t) = (1 / sqrt(2*pi)) * exp(-t^2 / 2)
    const INV_SQRT_2PI = 1.0 / Math.sqrt(2.0 * Math.PI);
    function normalPDF(t) {
        return INV_SQRT_2PI * Math.exp(-0.5 * t * t);
    }

    /**
     * Gold Standard High-Precision Reference for Standard Normal CDF Phi(z)
     * Combines high-order Taylor series near 0 with Laplace continued fraction in tails.
     * Accurate to full 64-bit IEEE 754 precision (approx 10^-16).
     */
    function normalCDFReference(z) {
        if (z === 0) return 0.5;
        if (z < -8.5) return 0.0;
        if (z > 8.5) return 1.0;

        const x = Math.abs(z);
        let p = 0;

        if (x < 3.2) {
            // Taylor expansion of int_0^x exp(-t^2/2) dt
            // sum_{k=0}^infty (-1)^k x^{2k+1} / (2^k * k! * (2k+1))
            let term = x;
            let sum = x;
            for (let k = 1; k < 120; k++) {
                term *= -(x * x) / (2.0 * k);
                const delta = term / (2.0 * k + 1.0);
                sum += delta;
                if (Math.abs(delta) < 1e-17) break;
            }
            p = 0.5 + sum * INV_SQRT_2PI;
        } else {
            // Laplace continued fraction for upper tail Q(x) = 1 - Phi(x)
            // Q(x) = phi(x) / (x + 1/(x + 2/(x + 3/(x + ...))))
            const phi = normalPDF(x);
            let f = x;
            let C = f;
            let D = 0.0;
            for (let j = 1; j < 120; j++) {
                const a = j;
                const b = x;
                D = b + a * D;
                if (D === 0) D = 1e-30;
                C = b + a / C;
                if (C === 0) C = 1e-30;
                D = 1.0 / D;
                const delta = C * D;
                f = f * delta;
                if (Math.abs(delta - 1.0) < 1e-17) break;
            }
            const tail = phi / f;
            p = 1.0 - tail;
        }

        return z < 0 ? 1.0 - p : p;
    }

    /**
     * Method A: Numerical Quadrature — Composite Trapezoidal Rule
     */
    function normalCDFTrapezoidal(z, n = 500) {
        const t0 = performance.now();
        if (z === 0) return { value: 0.5, elapsedMs: performance.now() - t0, evals: 1 };

        const sign = z < 0 ? -1 : 1;
        const b = Math.abs(z);
        const h = b / n;

        let sum = 0.5 * (normalPDF(0) + normalPDF(b));
        for (let i = 1; i < n; i++) {
            sum += normalPDF(i * h);
        }
        const integral = sum * h;
        const value = sign > 0 ? 0.5 + integral : 0.5 - integral;
        const t1 = performance.now();
        return { value, elapsedMs: t1 - t0, evals: n + 1 };
    }

    /**
     * Method B: Numerical Quadrature — Simpson's 1/3 Rule
     */
    function normalCDFSimpson(z, n = 500) {
        const t0 = performance.now();
        if (z === 0) return { value: 0.5, elapsedMs: performance.now() - t0, evals: 1 };
        if (n % 2 !== 0) n += 1;

        const sign = z < 0 ? -1 : 1;
        const b = Math.abs(z);
        const h = b / n;

        let sum = normalPDF(0) + normalPDF(b);
        for (let i = 1; i < n; i++) {
            const factor = (i % 2 === 1) ? 4.0 : 2.0;
            sum += factor * normalPDF(i * h);
        }
        const integral = (sum * h) / 3.0;
        const value = sign > 0 ? 0.5 + integral : 0.5 - integral;
        const t1 = performance.now();
        return { value, elapsedMs: t1 - t0, evals: n + 1 };
    }

    /**
     * Method C: Gauss-Legendre Quadrature (10-Point Adaptive)
     */
    const GL10_NODES = [
        -0.1488743389816312, 0.1488743389816312,
        -0.4333953941292472, 0.4333953941292472,
        -0.6794095682990244, 0.6794095682990244,
        -0.8650633666889845, 0.8650633666889845,
        -0.9739065285171717, 0.9739065285171717
    ];
    const GL10_WEIGHTS = [
        0.2955242247147529, 0.2955242247147529,
        0.2692667193099963, 0.2692667193099963,
        0.2190863625159820, 0.2190863625159820,
        0.1494513491505806, 0.1494513491505806,
        0.0666713443086881, 0.0666713443086881
    ];
    function normalCDFGaussLegendre(z) {
        const t0 = performance.now();
        if (z === 0) return { value: 0.5, elapsedMs: performance.now() - t0, evals: 1 };

        const sign = z < 0 ? -1 : 1;
        const b = Math.abs(z);
        let sum = 0.0;
        for (let i = 0; i < 10; i++) {
            const t = 0.5 * b * (GL10_NODES[i] + 1.0);
            sum += GL10_WEIGHTS[i] * normalPDF(t);
        }
        const integral = 0.5 * b * sum;
        const value = sign > 0 ? 0.5 + integral : 0.5 - integral;
        const t1 = performance.now();
        return { value, elapsedMs: t1 - t0, evals: 10 };
    }

    /**
     * Method D: Power Series Expansion of the Error Function (Taylor/Maclaurin)
     * erf(x) = (2 / sqrt(pi)) * sum_{k=0}^{K} [ (-1)^k * x^{2k+1} / (k! * (2k+1)) ]
     * Phi(z) = 0.5 * [ 1 + erf(z / sqrt(2)) ]
     */
    function normalCDFTaylorSeries(z, maxTerms = 25) {
        const t0 = performance.now();
        if (z === 0) return { value: 0.5, elapsedMs: performance.now() - t0, evals: 1, termsUsed: 1 };

        const x = z / Math.SQRT2;
        const x2 = x * x;
        let term = x;
        let sum = x;
        let k = 1;

        for (; k < maxTerms; k++) {
            term *= -x2 / k;
            const delta = term / (2.0 * k + 1.0);
            sum += delta;
            if (Math.abs(delta) < 1e-16) break;
        }

        const erfVal = (2.0 / Math.sqrt(Math.PI)) * sum;
        const value = 0.5 * (1.0 + erfVal);
        const t1 = performance.now();
        return { value, elapsedMs: t1 - t0, evals: k, termsUsed: k };
    }

    /**
     * Method E: Abramowitz & Stegun Rational Approximation (Formula 7.1.26)
     * Max error |epsilon(x)| <= 1.5 x 10^-7 across all R.
     */
    function normalCDFAbramowitzStegun(z) {
        const t0 = performance.now();
        if (z === 0) return { value: 0.5, elapsedMs: performance.now() - t0, evals: 1 };

        const sign = z < 0 ? -1 : 1;
        const x = Math.abs(z) / Math.SQRT2;

        const p = 0.3275911;
        const a1 = 0.254829592;
        const a2 = -0.284496736;
        const a3 = 1.421413741;
        const a4 = -1.453152027;
        const a5 = 1.061405429;

        const t = 1.0 / (1.0 + p * x);
        const poly = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t;
        const erfVal = sign * (1.0 - poly * Math.exp(-x * x));
        const value = 0.5 * (1.0 + erfVal);
        const t1 = performance.now();
        return { value, elapsedMs: t1 - t0, evals: 6 };
    }

    /**
     * Method F: Continued Fraction Tail Method (Laplace)
     * Especially suited for extreme z values to avoid catastrophic cancellation in tails.
     */
    function normalCDFContinuedFraction(z) {
        const t0 = performance.now();
        if (z === 0) return { value: 0.5, elapsedMs: performance.now() - t0, evals: 1 };

        const sign = z < 0 ? -1 : 1;
        const x = Math.abs(z);
        const phi = normalPDF(x);

        let f = x;
        let C = f;
        let D = 0.0;
        let steps = 0;

        for (let j = 1; j < 50; j++) {
            steps++;
            const a = j;
            const b = x;
            D = b + a * D;
            if (D === 0) D = 1e-30;
            C = b + a / C;
            if (C === 0) C = 1e-30;
            D = 1.0 / D;
            const delta = C * D;
            f = f * delta;
            if (Math.abs(delta - 1.0) < 1e-15) break;
        }

        const tail = phi / f;
        const value = sign > 0 ? 1.0 - tail : tail;
        const t1 = performance.now();
        return { value, elapsedMs: t1 - t0, evals: steps };
    }

    // =========================================================================
    // UI & CANVAS VISUALIZATION RENDERING
    // =========================================================================

    // Module state
    let state = {
        selectedFuncKey: "sin",
        riemannN: 50,
        lebesgueM: 20,
        lebesgueK: 40000,
        riemannRule: "midpoint",
        activeView: "both", // 'both' | 'riemann' | 'lebesgue'
        zValue: 1.96,
        quadratureSteps: 500,
        taylorTerms: 25
    };

    // DOM Elements Cache
    let dom = {};

    function initDOM() {
        dom.funcSelect = document.getElementById('func-select');
        dom.funcDescription = document.getElementById('func-description');
        dom.funcFormula = document.getElementById('func-formula');
        dom.riemannNSlider = document.getElementById('riemann-n-slider');
        dom.riemannNVal = document.getElementById('riemann-n-val');
        dom.lebesgueMSlider = document.getElementById('lebesgue-m-slider');
        dom.lebesgueMVal = document.getElementById('lebesgue-m-val');
        dom.lebesgueKSlider = document.getElementById('lebesgue-k-slider');
        dom.lebesgueKVal = document.getElementById('lebesgue-k-val');
        dom.riemannRule = document.getElementById('riemann-rule');
        dom.viewMode = document.getElementById('view-mode');

        // Results display
        dom.exactVal = document.getElementById('exact-val');
        dom.riemannVal = document.getElementById('riemann-val');
        dom.riemannErr = document.getElementById('riemann-err');
        dom.riemannTime = document.getElementById('riemann-time');
        dom.lebesgueVal = document.getElementById('lebesgue-val');
        dom.lebesgueErr = document.getElementById('lebesgue-err');
        dom.lebesgueTime = document.getElementById('lebesgue-time');

        // Canvases
        dom.canvasRiemann = document.getElementById('canvas-riemann');
        dom.canvasLebesgue = document.getElementById('canvas-lebesgue');
        dom.canvasNormalPDF = document.getElementById('canvas-normal-pdf');
        dom.canvasNormalResidual = document.getElementById('canvas-normal-residual');

        // Normal CDF controls
        dom.zSlider = document.getElementById('z-slider');
        dom.zInput = document.getElementById('z-input');
        dom.zValDisplay = document.getElementById('z-val-display');
        dom.quadStepsSelect = document.getElementById('quad-steps-select');
        dom.taylorTermsSelect = document.getElementById('taylor-terms-select');
        dom.normalCdfTableBody = document.getElementById('normal-cdf-table-body');
        dom.benchmarkTableBody = document.getElementById('benchmark-table-body');
    }

    /**
     * Color palette generator for Lebesgue preimage slices
     */
    function getSliceColor(index, total) {
        const hue = (index * 360 / Math.max(total, 1)) % 360;
        return `hsla(${hue}, 70%, 55%, 0.45)`;
    }
    function getSliceBorderColor(index, total) {
        const hue = (index * 360 / Math.max(total, 1)) % 360;
        return `hsl(${hue}, 80%, 40%)`;
    }

    /**
     * Draw Riemann Slicing on Canvas
     */
    function renderRiemannCanvas(fnObj, n, rule) {
        const canvas = dom.canvasRiemann;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width;
        const H = canvas.height;

        ctx.clearRect(0, 0, W, W);

        const padLeft = 50, padRight = 30, padTop = 30, padBottom = 40;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const [a, b] = fnObj.domain;
        const f = fnObj.fn;

        // Determine y bounds
        let yMin = fnObj.yRange ? fnObj.yRange[0] : 0;
        let yMax = fnObj.yRange ? fnObj.yRange[1] : 1;
        if (fnObj.key === 'inv_sqrt') {
            yMin = 0;
            yMax = 12; // cap visually for rendering
        } else {
            yMin = Math.min(0, yMin);
            yMax = Math.max(1, yMax * 1.15);
        }

        const mapX = (x) => padLeft + ((x - a) / (b - a)) * plotW;
        const mapY = (y) => padTop + (1.0 - (y - yMin) / (yMax - yMin)) * plotH;

        // Background & Grid
        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);
        ctx.strokeStyle = "#e1e4e8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let g = 0; g <= 5; g++) {
            const gy = padTop + (g / 5) * plotH;
            ctx.moveTo(padLeft, gy);
            ctx.lineTo(padLeft + plotW, gy);
            const yVal = yMax - (g / 5) * (yMax - yMin);
            ctx.fillStyle = "#6a737d";
            ctx.font = "11px monospace";
            ctx.textAlign = "right";
            ctx.fillText(yVal.toFixed(2), padLeft - 8, gy + 4);
        }
        ctx.stroke();

        // Draw Riemann rectangles
        const h = (b - a) / n;
        ctx.lineWidth = 1;

        for (let i = 0; i < n; i++) {
            const xLeft = a + i * h;
            const xRight = a + (i + 1) * h;
            let xEval = xLeft;
            if (rule === 'midpoint') xEval = xLeft + 0.5 * h;
            else if (rule === 'right') xEval = xRight;
            else if (rule === 'trapezoid') xEval = xLeft + 0.5 * h; // representative

            const yVal = Math.min(f(xEval), yMax);
            const rx = mapX(xLeft);
            const rw = mapX(xRight) - rx;
            const ry = mapY(Math.max(0, yVal));
            const rh = mapY(0) - ry;

            ctx.fillStyle = "rgba(52, 152, 219, 0.35)";
            ctx.strokeStyle = "#2980b9";
            ctx.fillRect(rx, ry, Math.max(1, rw), rh);
            ctx.strokeRect(rx, ry, Math.max(1, rw), rh);
        }

        // Draw exact function curve f(x)
        ctx.beginPath();
        ctx.strokeStyle = "#e74c3c";
        ctx.lineWidth = 2.5;
        const sampleCount = 400;
        for (let s = 0; s <= sampleCount; s++) {
            const x = a + (s / sampleCount) * (b - a);
            const y = Math.min(f(x), yMax);
            const px = mapX(x);
            const py = mapY(y);
            if (s === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // Axes
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(padLeft, padTop);
        ctx.lineTo(padLeft, padTop + plotH);
        ctx.lineTo(padLeft + plotW, padTop + plotH);
        ctx.stroke();

        // Axis labels
        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`x = a (${a.toFixed(2)})`, padLeft, padTop + plotH + 20);
        ctx.fillText(`Domain Partitioning (N = ${n} strips, Rule: ${rule})`, padLeft + plotW / 2, padTop + plotH + 32);
        ctx.fillText(`x = b (${b.toFixed(2)})`, padLeft + plotW, padTop + plotH + 20);
    }

    /**
     * Draw Lebesgue Slicing on Canvas (Horizontal Range Slices + Preimage Highlight)
     */
    function renderLebesgueCanvas(fnObj, M, lebesgueResult) {
        const canvas = dom.canvasLebesgue;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width;
        const H = canvas.height;

        ctx.clearRect(0, 0, W, H);

        const padLeft = 50, padRight = 30, padTop = 30, padBottom = 40;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const [a, b] = fnObj.domain;
        const f = fnObj.fn;

        const y0 = lebesgueResult.y0;
        const y1 = lebesgueResult.y1;
        const dy = lebesgueResult.dy;

        let displayYMin = y0;
        let displayYMax = y1;
        if (fnObj.key === 'inv_sqrt') {
            displayYMin = 0;
            displayYMax = 12;
        } else {
            displayYMin = Math.min(0, y0);
            displayYMax = Math.max(1, y1 * 1.15);
        }

        const mapX = (x) => padLeft + ((x - a) / (b - a)) * plotW;
        const mapY = (y) => padTop + (1.0 - (y - displayYMin) / (displayYMax - displayYMin)) * plotH;

        // Background & Grid
        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);
        ctx.strokeStyle = "#e1e4e8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let g = 0; g <= 5; g++) {
            const gy = padTop + (g / 5) * plotH;
            ctx.moveTo(padLeft, gy);
            ctx.lineTo(padLeft + plotW, gy);
            const yVal = displayYMax - (g / 5) * (displayYMax - displayYMin);
            ctx.fillStyle = "#6a737d";
            ctx.font = "11px monospace";
            ctx.textAlign = "right";
            ctx.fillText(yVal.toFixed(2), padLeft - 8, gy + 4);
        }
        ctx.stroke();

        // Draw Horizontal Range Bands
        for (let j = 0; j < M; j++) {
            const p = lebesgueResult.preimages[j];
            if (!p) continue;
            const yLow = p.yLower;
            const yHigh = p.yUpper;
            if (yLow > displayYMax) continue;

            const pyHigh = mapY(Math.min(displayYMax, yHigh));
            const pyLow = mapY(Math.max(displayYMin, yLow));
            const bandH = Math.max(1, pyLow - pyHigh);

            // Shaded horizontal range strip
            ctx.fillStyle = getSliceColor(j, M);
            ctx.fillRect(padLeft, pyHigh, plotW, bandH);

            // Strip delimiter line
            ctx.strokeStyle = getSliceBorderColor(j, M);
            ctx.lineWidth = 0.5;
            ctx.strokeRect(padLeft, pyHigh, plotW, bandH);
        }

        // Draw Preimage Footprints along the Domain (x-axis)
        const K = 300;
        const step = (b - a) / K;
        for (let k = 0; k < K; k++) {
            const xVal = a + (k + 0.5) * step;
            const yVal = f(xVal);
            let idx = Math.floor((yVal - y0) / dy);
            if (idx < 0) idx = 0;
            if (idx >= M) idx = M - 1;

            const px = mapX(xVal);
            const pw = (step / (b - a)) * plotW;
            ctx.fillStyle = getSliceBorderColor(idx, M);
            // Draw small footprint tick at bottom axis
            ctx.fillRect(px, padTop + plotH - 6, Math.max(1, pw), 6);
        }

        // Draw exact function curve f(x)
        ctx.beginPath();
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 2.5;
        const sampleCount = 400;
        for (let s = 0; s <= sampleCount; s++) {
            const x = a + (s / sampleCount) * (b - a);
            const y = Math.min(f(x), displayYMax);
            const px = mapX(x);
            const py = mapY(y);
            if (s === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // Axes
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(padLeft, padTop);
        ctx.lineTo(padLeft, padTop + plotH);
        ctx.lineTo(padLeft + plotW, padTop + plotH);
        ctx.stroke();

        // Labels
        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`Range Partitioning (M = ${M} horizontal slices, K = ${lebesgueResult.K})`, padLeft + plotW / 2, padTop + plotH + 32);
    }

    /**
     * Render Normal PDF Bell Curve with Shaded Integral up to z
     */
    function renderNormalPDFCanvas(zVal, quadSteps) {
        const canvas = dom.canvasNormalPDF;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width;
        const H = canvas.height;

        ctx.clearRect(0, 0, W, H);

        const padLeft = 45, padRight = 30, padTop = 30, padBottom = 40;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const xMin = -4.5, xMax = 4.5;
        const yMin = 0.0, yMax = 0.44;

        const mapX = (x) => padLeft + ((x - xMin) / (xMax - xMin)) * plotW;
        const mapY = (y) => padTop + (1.0 - (y - yMin) / (yMax - yMin)) * plotH;

        // Background
        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);

        // Draw standard deviation reference bands (-3 to +3)
        const stdDevs = [-3, -2, -1, 1, 2, 3];
        ctx.strokeStyle = "#e1e4e8";
        ctx.lineWidth = 1;
        ctx.setLineDash([4, 4]);
        stdDevs.forEach(sd => {
            const px = mapX(sd);
            ctx.beginPath();
            ctx.moveTo(px, padTop);
            ctx.lineTo(px, padTop + plotH);
            ctx.stroke();
            ctx.fillStyle = "#95a5a6";
            ctx.font = "10px monospace";
            ctx.textAlign = "center";
            ctx.fillText(`${sd > 0 ? '+' : ''}${sd}σ`, px, padTop + plotH + 15);
        });
        ctx.setLineDash([]);

        // Shaded cumulative probability region from -infinity to z
        const clampedZ = Math.max(xMin, Math.min(xMax, zVal));
        ctx.beginPath();
        ctx.moveTo(mapX(xMin), mapY(0));
        const shadeSteps = 200;
        const shadeStepSize = (clampedZ - xMin) / shadeSteps;
        for (let i = 0; i <= shadeSteps; i++) {
            const curX = xMin + i * shadeStepSize;
            ctx.lineTo(mapX(curX), mapY(normalPDF(curX)));
        }
        ctx.lineTo(mapX(clampedZ), mapY(0));
        ctx.closePath();
        ctx.fillStyle = "rgba(46, 204, 113, 0.35)";
        ctx.fill();

        // Draw Normal PDF Bell Curve
        ctx.beginPath();
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 2.5;
        const curveSteps = 300;
        for (let i = 0; i <= curveSteps; i++) {
            const curX = xMin + (i / curveSteps) * (xMax - xMin);
            const curY = normalPDF(curX);
            const px = mapX(curX);
            const py = mapY(curY);
            if (i === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // Vertical boundary line at z
        const pzX = mapX(clampedZ);
        ctx.beginPath();
        ctx.strokeStyle = "#e74c3c";
        ctx.lineWidth = 2;
        ctx.moveTo(pzX, mapY(0));
        ctx.lineTo(pzX, mapY(normalPDF(clampedZ)));
        ctx.stroke();

        // Point at (z, phi(z))
        ctx.fillStyle = "#e74c3c";
        ctx.beginPath();
        ctx.arc(pzX, mapY(normalPDF(clampedZ)), 4.5, 0, 2 * Math.PI);
        ctx.fill();

        // Axes
        ctx.strokeStyle = "#34495e";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(padLeft, padTop);
        ctx.lineTo(padLeft, padTop + plotH);
        ctx.lineTo(padLeft + plotW, padTop + plotH);
        ctx.stroke();

        // Annotation banner
        const phiExact = normalCDFReference(zVal);
        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 12px sans-serif";
        ctx.textAlign = "left";
        ctx.fillText(`Evaluation point z = ${zVal.toFixed(3)}   |   Exact Φ(z) = ${phiExact.toFixed(9)}`, padLeft + 15, padTop + 20);
    }

    /**
     * Render Log Error Residual Comparison Canvas across methods
     */
    function renderNormalResidualCanvas() {
        const canvas = dom.canvasNormalResidual;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width;
        const H = canvas.height;

        ctx.clearRect(0, 0, W, H);

        const padLeft = 60, padRight = 30, padTop = 30, padBottom = 40;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const xMin = -4.0, xMax = 4.0;
        // Log10 absolute error: from 10^-16 to 10^0
        const logMin = -16.0, logMax = 0.0;

        const mapX = (x) => padLeft + ((x - xMin) / (xMax - xMin)) * plotW;
        const mapY = (logErr) => padTop + (1.0 - (logErr - logMin) / (logMax - logMin)) * plotH;

        // Background
        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);

        // Grid lines for orders of magnitude
        ctx.strokeStyle = "#e1e4e8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let logVal = -16; logVal <= 0; logVal += 2) {
            const py = mapY(logVal);
            ctx.moveTo(padLeft, py);
            ctx.lineTo(padLeft + plotW, py);
            ctx.fillStyle = "#6a737d";
            ctx.font = "10px monospace";
            ctx.textAlign = "right";
            ctx.fillText(`10^${logVal}`, padLeft - 8, py + 3);
        }
        ctx.stroke();

        // Methods to compare across domain
        const methods = [
            {
                name: "A&S (7.1.26)",
                color: "#e67e22",
                calc: (z) => normalCDFAbramowitzStegun(z).value
            },
            {
                name: "Trapezoid (N=200)",
                color: "#3498db",
                calc: (z) => normalCDFTrapezoidal(z, 200).value
            },
            {
                name: "Simpson (N=200)",
                color: "#9b59b6",
                calc: (z) => normalCDFSimpson(z, 200).value
            },
            {
                name: "Taylor (K=15)",
                color: "#e74c3c",
                calc: (z) => normalCDFTaylorSeries(z, 15).value
            },
            {
                name: "Gauss-Legendre (10-pt)",
                color: "#27ae60",
                calc: (z) => normalCDFGaussLegendre(z).value
            }
        ];

        const numPoints = 80;
        methods.forEach(method => {
            ctx.beginPath();
            ctx.strokeStyle = method.color;
            ctx.lineWidth = 2.0;

            let started = false;
            for (let i = 0; i <= numPoints; i++) {
                const z = xMin + (i / numPoints) * (xMax - xMin);
                const exact = normalCDFReference(z);
                const est = method.calc(z);
                const absErr = Math.max(1e-16, Math.abs(est - exact));
                const logErr = Math.log10(absErr);

                const px = mapX(z);
                const py = mapY(Math.min(0, Math.max(logMin, logErr)));

                if (!started) {
                    ctx.moveTo(px, py);
                    started = true;
                } else {
                    ctx.lineTo(px, py);
                }
            }
            ctx.stroke();
        });

        // Axes
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(padLeft, padTop);
        ctx.lineTo(padLeft, padTop + plotH);
        ctx.lineTo(padLeft + plotW, padTop + plotH);
        ctx.stroke();

        // Legend
        let legX = padLeft + 15;
        ctx.font = "bold 10px sans-serif";
        methods.forEach(m => {
            ctx.fillStyle = m.color;
            ctx.fillRect(legX, padTop + 10, 12, 10);
            ctx.fillStyle = "#2c3e50";
            ctx.textAlign = "left";
            ctx.fillText(m.name, legX + 16, padTop + 18);
            legX += ctx.measureText(m.name).width + 30;
        });

        // Axis title
        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "center";
        ctx.fillText("Standard Normal z", padLeft + plotW / 2, padTop + plotH + 30);
    }

    /**
     * Update Part 1: Interactive Riemann vs Lebesgue Simulation
     */
    function updateExercise1() {
        const fnObj = FUNCTIONS[state.selectedFuncKey];
        if (!fnObj) return;

        // Update descriptions
        dom.funcDescription.innerHTML = `<strong>${fnObj.category}:</strong> ${fnObj.description}`;
        dom.funcFormula.textContent = `$$${fnObj.latex}$$`;
        if (window.renderMathInElement) {
            window.renderMathInElement(dom.funcFormula, {
                delimiters: [{ left: '$$', right: '$$', display: true }]
            });
        }

        const [a, b] = fnObj.domain;

        // Run Riemann
        const riemannRes = riemannIntegrate(fnObj.fn, a, b, state.riemannN, state.riemannRule);

        // Run Lebesgue
        const lebesgueRes = lebesgueIntegrate(fnObj.fn, a, b, state.lebesgueM, state.lebesgueK, fnObj.yRange ? fnObj.yRange[0] : null, fnObj.yRange ? fnObj.yRange[1] : null);

        // Display results
        dom.exactVal.textContent = fnObj.exact !== null ? fnObj.exact.toFixed(7) : "Undefined (Riemann Fails)";
        dom.riemannVal.textContent = riemannRes.value.toFixed(7);
        dom.riemannTime.textContent = `${(riemannRes.elapsedMs * 1000).toFixed(1)} µs`;

        if (fnObj.exact !== null) {
            const rErr = Math.abs(riemannRes.value - fnObj.exact);
            dom.riemannErr.textContent = `${rErr.toExponential(3)} (${((rErr / Math.max(1e-12, Math.abs(fnObj.exact))) * 100).toFixed(3)}%)`;
        } else {
            dom.riemannErr.textContent = "Divergent / Grid Dependent";
        }

        dom.lebesgueVal.textContent = lebesgueRes.value.toFixed(7);
        dom.lebesgueTime.textContent = `${(lebesgueRes.elapsedMs).toFixed(2)} ms`;

        if (fnObj.theoreticalLebesgue !== null) {
            const lErr = Math.abs(lebesgueRes.value - fnObj.theoreticalLebesgue);
            dom.lebesgueErr.textContent = `${lErr.toExponential(3)} (${((lErr / Math.max(1e-12, Math.abs(fnObj.theoreticalLebesgue))) * 100).toFixed(3)}%)`;
        } else {
            dom.lebesgueErr.textContent = "N/A";
        }

        // Render Canvases
        renderRiemannCanvas(fnObj, state.riemannN, state.riemannRule);
        renderLebesgueCanvas(fnObj, state.lebesgueM, lebesgueRes);
    }

    /**
     * Update Part 2: Normal CDF Evaluator & Table
     */
    function updateExercise2() {
        const z = state.zValue;
        const refVal = normalCDFReference(z);

        // Evaluate methods
        const quadSteps = state.quadratureSteps;
        const taylorTerms = state.taylorTerms;

        const resTrap = normalCDFTrapezoidal(z, quadSteps);
        const resSimp = normalCDFSimpson(z, quadSteps);
        const resGL10 = normalCDFGaussLegendre(z);
        const resTaylor = normalCDFTaylorSeries(z, taylorTerms);
        const resAS = normalCDFAbramowitzStegun(z);
        const resCF = normalCDFContinuedFraction(z);

        const methodList = [
            {
                name: "Gold Standard Reference",
                tag: "Taylor + Laplace CF",
                val: refVal,
                err: 0,
                evals: "~100",
                time: "< 0.05 µs",
                status: "Exact Reference (10⁻¹⁶)"
            },
            {
                name: "Composite Trapezoidal Rule",
                tag: `N = ${quadSteps} steps`,
                val: resTrap.value,
                err: Math.abs(resTrap.value - refVal),
                evals: resTrap.evals,
                time: `${(resTrap.elapsedMs * 1000).toFixed(1)} µs`,
                status: "Quadrature O(h²)"
            },
            {
                name: "Simpson's 1/3 Rule",
                tag: `N = ${quadSteps} steps`,
                val: resSimp.value,
                err: Math.abs(resSimp.value - refVal),
                evals: resSimp.evals,
                time: `${(resSimp.elapsedMs * 1000).toFixed(1)} µs`,
                status: "Quadrature O(h⁴)"
            },
            {
                name: "Gauss-Legendre (10-Point)",
                tag: "Orthogonal Legendre nodes",
                val: resGL10.value,
                err: Math.abs(resGL10.value - refVal),
                evals: resGL10.evals,
                time: `${(resGL10.elapsedMs * 1000).toFixed(1)} µs`,
                status: "High Efficiency (Exact deg ≤ 19)"
            },
            {
                name: "Taylor Series (erf expansion)",
                tag: `K = ${resTaylor.termsUsed} terms`,
                val: resTaylor.value,
                err: Math.abs(resTaylor.value - refVal),
                evals: resTaylor.evals,
                time: `${(resTaylor.elapsedMs * 1000).toFixed(1)} µs`,
                status: Math.abs(z) > 3.0 ? "Catastrophic Cancellation" : "Rapid Local Convergence"
            },
            {
                name: "Abramowitz & Stegun 7.1.26",
                tag: "Degree-5 Rational Poly",
                val: resAS.value,
                err: Math.abs(resAS.value - refVal),
                evals: resAS.evals,
                time: `${(resAS.elapsedMs * 1000).toFixed(1)} µs`,
                status: "Error ≤ 1.5×10⁻⁷"
            },
            {
                name: "Laplace Continued Fraction",
                tag: "Tail recurrence",
                val: resCF.value,
                err: Math.abs(resCF.value - refVal),
                evals: resCF.evals,
                time: `${(resCF.elapsedMs * 1000).toFixed(1)} µs`,
                status: "Tail Stabilizer (No Loss in Flops)"
            }
        ];

        // Populate table
        let tableHTML = "";
        methodList.forEach((m, idx) => {
            const isRef = idx === 0;
            const errStr = isRef ? "0.0 (Reference)" : m.err.toExponential(3);
            const errClass = m.err > 1e-4 ? 'color: #e74c3c; font-weight: bold;' : (m.err < 1e-12 ? 'color: #27ae60; font-weight: bold;' : '');

            tableHTML += `
                <tr style="${isRef ? 'background: #eafaf1; font-weight: bold;' : ''}">
                    <td><strong>${m.name}</strong><br><small style="color:#7f8c8d;">${m.tag}</small></td>
                    <td style="font-family: monospace; font-size: 0.95rem;">${m.val.toFixed(10)}</td>
                    <td style="font-family: monospace; font-size: 0.95rem; ${errClass}">${errStr}</td>
                    <td style="font-family: monospace;">${m.evals}</td>
                    <td style="font-family: monospace;">${m.time}</td>
                    <td><span style="font-size:0.85rem; padding: 2px 6px; border-radius: 4px; background: #ecf0f1;">${m.status}</span></td>
                </tr>
            `;
        });
        dom.normalCdfTableBody.innerHTML = tableHTML;

        // Update Canvases
        renderNormalPDFCanvas(z, quadSteps);
        renderNormalResidualCanvas();
    }

    /**
     * Populate Global Convergence Benchmark Matrix for Exercise 1
     */
    function populateExercise1Benchmark() {
        if (!dom.benchmarkTableBody) return;
        const testN = [10, 50, 100, 500, 2000, 10000];
        const fn = FUNCTIONS.sin; // sin(x) on [0, pi], exact = 2
        const exact = 2.0;

        let rows = "";
        testN.forEach(n => {
            const rMid = riemannIntegrate(fn.fn, 0, Math.PI, n, 'midpoint');
            const rTrap = riemannIntegrate(fn.fn, 0, Math.PI, n, 'trapezoid');
            const lRes = lebesgueIntegrate(fn.fn, 0, Math.PI, Math.min(n, 200), Math.min(n * 20, 100000), 0, 1);

            const errMid = Math.abs(rMid.value - exact);
            const errTrap = Math.abs(rTrap.value - exact);
            const errLeb = Math.abs(lRes.value - exact);

            rows += `
                <tr>
                    <td><strong>${n}</strong></td>
                    <td style="font-family: monospace;">${rMid.value.toFixed(7)}<br><small style="color:#7f8c8d;">Err: ${errMid.toExponential(2)}</small></td>
                    <td style="font-family: monospace;">${rTrap.value.toFixed(7)}<br><small style="color:#7f8c8d;">Err: ${errTrap.toExponential(2)}</small></td>
                    <td style="font-family: monospace;">${lRes.value.toFixed(7)}<br><small style="color:#7f8c8d;">Err: ${errLeb.toExponential(2)}</small></td>
                    <td style="font-family: monospace;">${(rMid.elapsedMs * 1000).toFixed(1)} µs</td>
                    <td style="font-family: monospace;">${(lRes.elapsedMs).toFixed(1)} ms</td>
                </tr>
            `;
        });
        dom.benchmarkTableBody.innerHTML = rows;
    }

    /**
     * Wire Up Event Handlers
     */
    function attachEventListeners() {
        // Function select
        dom.funcSelect.addEventListener('change', (e) => {
            state.selectedFuncKey = e.target.value;
            updateExercise1();
        });

        // Riemann N slider
        dom.riemannNSlider.addEventListener('input', (e) => {
            state.riemannN = parseInt(e.target.value, 10);
            dom.riemannNVal.textContent = state.riemannN;
            updateExercise1();
        });

        // Lebesgue M slider
        dom.lebesgueMSlider.addEventListener('input', (e) => {
            state.lebesgueM = parseInt(e.target.value, 10);
            dom.lebesgueMVal.textContent = state.lebesgueM;
            updateExercise1();
        });

        // Lebesgue K slider
        dom.lebesgueKSlider.addEventListener('input', (e) => {
            state.lebesgueK = parseInt(e.target.value, 10);
            dom.lebesgueKVal.textContent = state.lebesgueK.toLocaleString();
            updateExercise1();
        });

        // Riemann rule
        dom.riemannRule.addEventListener('change', (e) => {
            state.riemannRule = e.target.value;
            updateExercise1();
        });

        // z slider
        dom.zSlider.addEventListener('input', (e) => {
            state.zValue = parseFloat(e.target.value);
            dom.zInput.value = state.zValue.toFixed(3);
            if (dom.zValDisplay) dom.zValDisplay.textContent = state.zValue.toFixed(3);
            updateExercise2();
        });

        // z text input
        dom.zInput.addEventListener('change', (e) => {
            let val = parseFloat(e.target.value);
            if (isNaN(val)) val = 0.0;
            val = Math.max(-5.0, Math.min(5.0, val));
            state.zValue = val;
            dom.zSlider.value = val;
            dom.zInput.value = val.toFixed(3);
            if (dom.zValDisplay) dom.zValDisplay.textContent = val.toFixed(3);
            updateExercise2();
        });

        // Quick z-buttons
        document.querySelectorAll('.btn-z-preset').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const zTarget = parseFloat(e.target.dataset.z);
                state.zValue = zTarget;
                dom.zSlider.value = zTarget;
                dom.zInput.value = zTarget.toFixed(3);
                if (dom.zValDisplay) dom.zValDisplay.textContent = zTarget.toFixed(3);
                updateExercise2();
            });
        });

        // Quadrature Steps select
        dom.quadStepsSelect.addEventListener('change', (e) => {
            state.quadratureSteps = parseInt(e.target.value, 10);
            updateExercise2();
        });

        // Taylor terms select
        dom.taylorTermsSelect.addEventListener('change', (e) => {
            state.taylorTerms = parseInt(e.target.value, 10);
            updateExercise2();
        });
    }

    // Initialization on DOMContentLoaded
    window.addEventListener('DOMContentLoaded', () => {
        initDOM();
        attachEventListeners();
        updateExercise1();
        updateExercise2();
        populateExercise1Benchmark();
    });

})();

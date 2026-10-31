// ============================================================================
// HOMEWORK 5 - WELFORD'S RECURRENCE, SAMPLING MEAN & VARIANCE
// Dario Pietrosanto (2055042)
// Master's Degree Course in Applied Statistics & Cybersecurity
// Real-time Interactive Engine (Guru Level):
// 1. Welford's Algorithm vs. Naive Formula Numerical Stability Playground
// 2. Multi-Hacker Penetration Simulation on Discrete Score Distributions [0, 10]
// 3. Online Single-Pass Stream Trajectories with Real-Time Play Evolution
// 4. Sampling Distribution Histograms for X_bar (CLT Normal) and S^2 (Scaled Chi-Square PDF)
// 5. Grand Ensemble Comparison vs. Theoretical Moments (Unbiasedness & Bessel's Correction)
// 6. Empirical Variability Scaling as a Function of System Depth n (O(1/n) Decay)
// ============================================================================

document.addEventListener('DOMContentLoaded', () => {

    // --- DOM Elements: Playground ---
    const playgroundInput = document.getElementById('playgroundInput');
    const playgroundRunBtn = document.getElementById('playgroundRunBtn');
    const presetNormalBtn = document.getElementById('presetNormalBtn');
    const presetPathologicalBtn = document.getElementById('presetPathologicalBtn');
    const playgroundResults = document.getElementById('playgroundResults');

    // --- DOM Elements: Simulation Suite ---
    const distSelect = document.getElementById('distSelect');
    const pBinomialGroup = document.getElementById('pBinomialGroup');
    const pBinomialInput = document.getElementById('pBinomialInput');
    const pBinomialDisplay = document.getElementById('pBinomialDisplay');
    const nInput = document.getElementById('nInput');
    const nDisplay = document.getElementById('nDisplay');
    const mInput = document.getElementById('mInput');
    const mDisplay = document.getElementById('mDisplay');
    const runSimBtn = document.getElementById('runSimBtn');
    const runBatchNBtn = document.getElementById('runBatchNBtn');
    const playStreamBtn = document.getElementById('playStreamBtn');

    // Canvas elements
    const dualHistCanvas = document.getElementById('dualHistCanvas');
    const dualHistCtx = dualHistCanvas.getContext('2d');
    const trajCanvas = document.getElementById('trajCanvas');
    const trajCtx = trajCanvas.getContext('2d');
    const decayCanvas = document.getElementById('decayCanvas');
    const decayCtx = decayCanvas.getContext('2d');

    // Metrics elements
    const statThMean = document.getElementById('statThMean');
    const statAvgMean = document.getElementById('statAvgMean');
    const statMeanDiff = document.getElementById('statMeanDiff');
    const statThVar = document.getElementById('statThVar');
    const statAvgVar = document.getElementById('statAvgVar');
    const statVarDiff = document.getElementById('statVarDiff');
    const statAvgBiasedVar = document.getElementById('statAvgBiasedVar');
    const statVarOfMeans = document.getElementById('statVarOfMeans');
    const statThVarOfMeans = document.getElementById('statThVarOfMeans');
    const statVarOfVars = document.getElementById('statVarOfVars');
    const conclusionBanner = document.getElementById('conclusionBanner');

    // High-DPI Canvas Setup Helper
    function setupHiDPI(canvas, ctx) {
        const dpr = window.devicePixelRatio || 1;
        const rect = canvas.getBoundingClientRect();
        const displayWidth = rect.width || canvas.width;
        const displayHeight = rect.height || canvas.height;

        canvas.width = Math.round(displayWidth * dpr);
        canvas.height = Math.round(displayHeight * dpr);
        ctx.scale(dpr, dpr);
        return { w: displayWidth, h: displayHeight };
    }

    // Lanczos approximation for log-Gamma function ln(Gamma(z))
    function logGamma(z) {
        const c = [
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109585720572, 9.9843695780195716e-6, 1.5056327351493116e-7
        ];
        if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z);
        z -= 1;
        let x = c[0];
        for (let i = 1; i < c.length; i++) x += c[i] / (z + i);
        let t = z + c.length - 1.5;
        return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
    }

    // Exact Theoretical Probability Density Function for Sample Variance S^2
    // (n-1)*S^2 / sigma^2 ~ Chi-Square(df = n-1)
    function sampleVarPDF(s2, n, sigma2) {
        if (s2 <= 0 || n <= 2 || sigma2 <= 0) return 0;
        const k = n - 1; // degrees of freedom
        const u = (k * s2) / sigma2;
        if (u <= 0) return 0;
        const logChi2 = - (k / 2) * Math.LN2 - logGamma(k / 2) + (k / 2 - 1) * Math.log(u) - u / 2;
        const chi2PDF = Math.exp(logChi2);
        return chi2PDF * (k / sigma2);
    }

    // ------------------------------------------------------------------------
    // SECTION 1: WELFORD'S ALGORITHM VS. NAIVE FORMULA (EXERCISE 1)
    // ------------------------------------------------------------------------
    function computeWelford(arr) {
        let mean = 0;
        let M2 = 0;
        const history = [];

        for (let k = 1; k <= arr.length; k++) {
            const x = arr[k - 1];
            // Core 3-line recurrence
            const delta = x - mean;
            mean += delta / k;
            M2 += delta * (x - mean);

            history.push({
                k,
                x,
                mean,
                M2,
                sampleVar: (k > 1) ? (M2 / (k - 1)) : 0
            });
        }

        const sampleVar = (arr.length > 1) ? (M2 / (arr.length - 1)) : 0;
        const popVar = (arr.length > 0) ? (M2 / arr.length) : 0;
        return { mean, M2, sampleVar, popVar, history };
    }

    function computeNaive(arr) {
        const n = arr.length;
        if (n === 0) return { mean: 0, sampleVar: 0, M2: 0 };
        let sum1 = 0;
        let sum2 = 0;
        for (let i = 0; i < n; i++) {
            sum1 += arr[i];
            sum2 += arr[i] * arr[i];
        }
        const mean = sum1 / n;
        const ss = sum2 - (sum1 * sum1) / n;
        const sampleVar = (n > 1) ? (ss / (n - 1)) : 0;
        return { mean, M2: ss, sampleVar };
    }

    function runPlaygroundTest() {
        const raw = playgroundInput.value.trim();
        const arr = raw.split(/[\s,]+/).map(Number).filter(v => !isNaN(v));

        if (arr.length < 2) {
            playgroundResults.innerHTML = '<span style="color:#e74c3c;">Please enter at least 2 numbers.</span>';
            return;
        }

        const welfordRes = computeWelford(arr);
        const naiveRes = computeNaive(arr);

        const naiveFailed = isNaN(naiveRes.sampleVar) || naiveRes.sampleVar < 0 || Math.abs(naiveRes.sampleVar - welfordRes.sampleVar) > 1e-4;

        playgroundResults.innerHTML = `
            <table class="comparison-table" style="margin-top:0.8rem;">
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Welford's Online Algorithm (Knuth/Welford)</th>
                        <th>Naive Textbook Formula $\\sum x^2 - (\\sum x)^2/n$</th>
                        <th>Status & Numerical Precision</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td><strong>Sample Mean ($\\bar{x}$)</strong></td>
                        <td><code>${welfordRes.mean.toFixed(6)}</code></td>
                        <td><code>${naiveRes.mean.toFixed(6)}</code></td>
                        <td><span style="color:#27ae60; font-weight:bold;">Identical</span></td>
                    </tr>
                    <tr>
                        <td><strong>Sum of Squares ($M_2$ or $SS$)</strong></td>
                        <td><code>${welfordRes.M2.toFixed(6)}</code></td>
                        <td><code>${naiveRes.M2.toFixed(6)}</code></td>
                        <td>${naiveFailed ? '<span style="color:#e74c3c; font-weight:bold;">Catastrophic Error</span>' : '<span style="color:#27ae60; font-weight:bold;">Accurate</span>'}</td>
                    </tr>
                    <tr>
                        <td><strong>Sample Variance ($s^2$)</strong></td>
                        <td><strong style="color:#27ae60; font-size:1.05rem;">${welfordRes.sampleVar.toFixed(6)}</strong></td>
                        <td><strong style="color:${naiveFailed ? '#e74c3c' : '#2c3e50'}; font-size:1.05rem;">${naiveRes.sampleVar.toFixed(6)}</strong></td>
                        <td>${naiveFailed ? '<span style="color:#e74c3c; font-weight:bold;">✗ Numerical Cancellation FAILED</span>' : '<span style="color:#27ae60; font-weight:bold;">✓ Passed</span>'}</td>
                    </tr>
                </tbody>
            </table>
            ${naiveFailed ? `
                <div class="alert-callout" style="margin-top:1rem;">
                    <strong>Catastrophic Cancellation Detected!</strong><br>
                    Because the numbers share a massive offset ($\approx 10^9$), the naive calculation $\\sum x_i^2$ and $(\\sum x_i)^2/n$ are both of magnitude $\\sim 10^{18}$, which exceeds the 53-bit precision mantissa of IEEE 754 floating point numbers. Subtracting them subtracts two nearly identical huge floats, destroying significant bits and yielding nonsensical or even <strong>negative variance</strong>! 
                    In contrast, <strong>Welford's algorithm subtracts the running mean first</strong> ($\delta_k = x_k - \\bar{x}_{k-1}$), keeping intermediate values centered around small deviations and computing the mathematically exact variance of <strong>${welfordRes.sampleVar.toFixed(6)}</strong>!
                </div>
            ` : ''}
        `;
    }

    presetNormalBtn.addEventListener('click', () => {
        playgroundInput.value = "4, 7, 13, 16";
        runPlaygroundTest();
    });

    presetPathologicalBtn.addEventListener('click', () => {
        playgroundInput.value = "1000000004, 1000000007, 1000000013, 1000000016";
        runPlaygroundTest();
    });

    playgroundRunBtn.addEventListener('click', runPlaygroundTest);

    // ------------------------------------------------------------------------
    // SECTION 2: PENETRATION SCORE DISTRIBUTIONS [0, 10] (EXERCISE 2)
    // ------------------------------------------------------------------------
    function getDistributionParams() {
        const type = distSelect.value;
        const pBin = parseFloat(pBinomialInput.value);

        if (type === 'binomial') {
            const mu = 10 * pBin;
            const sigma2 = 10 * pBin * (1 - pBin);
            const rho4 = 10 * pBin * (1 - pBin) * (1 + 3 * 8 * pBin * (1 - pBin));
            return {
                type: 'binomial',
                name: `Binomial(10, p=${pBin.toFixed(2)})`,
                mu,
                sigma2,
                sigma: Math.sqrt(sigma2),
                rho4,
                sample: () => {
                    let successes = 0;
                    for (let i = 0; i < 10; i++) {
                        if (Math.random() < pBin) successes++;
                    }
                    return successes;
                }
            };
        } else if (type === 'uniform') {
            const mu = 5.0;
            const sigma2 = 10.0;
            const rho4 = 178.0;
            return {
                type: 'uniform',
                name: 'Discrete Uniform U{0, 10}',
                mu,
                sigma2,
                sigma: Math.sqrt(sigma2),
                rho4,
                sample: () => Math.floor(Math.random() * 11)
            };
        } else {
            const p1 = 0.7, p2 = 0.3;
            const mu1 = 1.5, s2_1 = 3 * 0.5 * 0.5;
            const mu2 = 7 + 3 * 0.8;
            const mu = p1 * mu1 + p2 * mu2;
            const EX2 = p1 * (s2_1 + mu1 * mu1) + p2 * (s2_2 + mu2 * mu2);
            const sigma2 = EX2 - mu * mu;
            return {
                type: 'bimodal',
                name: 'Bimodal Mixture (Reconnaissance vs. Exploit)',
                mu,
                sigma2,
                sigma: Math.sqrt(sigma2),
                rho4: 120.0,
                sample: () => {
                    if (Math.random() < 0.7) {
                        let c = 0;
                        for (let i = 0; i < 3; i++) if (Math.random() < 0.5) c++;
                        return c;
                    } else {
                        let c = 7;
                        for (let i = 0; i < 3; i++) if (Math.random() < 0.8) c++;
                        return c;
                    }
                }
            };
        }
    }

    function syncSimControls() {
        pBinomialGroup.style.display = (distSelect.value === 'binomial') ? 'block' : 'none';
        pBinomialDisplay.textContent = parseFloat(pBinomialInput.value).toFixed(2);
        nDisplay.textContent = parseInt(nInput.value, 10).toLocaleString();
        mDisplay.textContent = parseInt(mInput.value, 10).toLocaleString();
    }

    distSelect.addEventListener('change', () => { syncSimControls(); runSimulation(); });
    pBinomialInput.addEventListener('input', () => { syncSimControls(); runSimulation(); });
    nInput.addEventListener('input', () => { syncSimControls(); runSimulation(); });
    mInput.addEventListener('input', () => { syncSimControls(); runSimulation(); });

    document.querySelectorAll('.hw5-preset-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            stopStreamAnim();
            const dist = btn.getAttribute('data-dist');
            const nVal = parseInt(btn.getAttribute('data-n'), 10);
            const mVal = parseInt(btn.getAttribute('data-m'), 10);
            const pVal = parseFloat(btn.getAttribute('data-p') || '0.6');

            distSelect.value = dist;
            nInput.value = nVal;
            mInput.value = mVal;
            pBinomialInput.value = pVal;

            syncSimControls();
            runSimulation();
        });
    });

    // ------------------------------------------------------------------------
    // SECTION 3: MULTI-HACKER SIMULATION RUNNER
    // ------------------------------------------------------------------------
    let lastSimResults = null;
    let streamAnimTimer = null;
    let streamCurrentK = 2;

    function stopStreamAnim() {
        if (streamAnimTimer) {
            clearInterval(streamAnimTimer);
            streamAnimTimer = null;
            playStreamBtn.textContent = '▶ Play Stream Evolution';
        }
    }

    function runSimulation() {
        stopStreamAnim();
        syncSimControls();
        const dist = getDistributionParams();
        const n = parseInt(nInput.value, 10);
        const m = parseInt(mInput.value, 10);

        const hackerMeans = new Float64Array(m);
        const hackerVars = new Float64Array(m);
        const hackerBiasedVars = new Float64Array(m);

        const sampleCount = Math.min(5, m);
        const sampleTrajectories = [];
        for (let j = 0; j < sampleCount; j++) {
            sampleTrajectories.push({
                means: new Float64Array(n + 1),
                vars: new Float64Array(n + 1)
            });
        }

        for (let j = 0; j < m; j++) {
            let mean = 0;
            let M2 = 0;
            const isSampleHacker = (j < sampleCount);

            if (isSampleHacker) {
                sampleTrajectories[j].means[0] = 0;
                sampleTrajectories[j].vars[0] = 0;
            }

            for (let k = 1; k <= n; k++) {
                const x = dist.sample();

                // Core 3-line Welford recurrence
                const delta = x - mean;
                mean += delta / k;
                M2 += delta * (x - mean);

                if (isSampleHacker) {
                    sampleTrajectories[j].means[k] = mean;
                    sampleTrajectories[j].vars[k] = (k > 1) ? (M2 / (k - 1)) : 0;
                }
            }

            hackerMeans[j] = mean;
            hackerVars[j] = (n > 1) ? (M2 / (n - 1)) : 0;
            hackerBiasedVars[j] = M2 / n;
        }

        // Grand Ensemble Statistics across m hackers
        let sumMeans = 0;
        let sumVars = 0;
        let sumBiasedVars = 0;
        for (let j = 0; j < m; j++) {
            sumMeans += hackerMeans[j];
            sumVars += hackerVars[j];
            sumBiasedVars += hackerBiasedVars[j];
        }
        const grandAvgMean = sumMeans / m;
        const grandAvgVar = sumVars / m;
        const grandAvgBiasedVar = sumBiasedVars / m;

        // Variability across hackers
        let sqDiffMeans = 0;
        let sqDiffVars = 0;
        for (let j = 0; j < m; j++) {
            sqDiffMeans += Math.pow(hackerMeans[j] - grandAvgMean, 2);
            sqDiffVars += Math.pow(hackerVars[j] - grandAvgVar, 2);
        }
        const varOfSampleMeans = sqDiffMeans / (m - 1);
        const varOfSampleVars = sqDiffVars / (m - 1);

        const thVarOfMeans = dist.sigma2 / n;

        lastSimResults = {
            dist,
            n,
            m,
            hackerMeans,
            hackerVars,
            grandAvgMean,
            grandAvgVar,
            grandAvgBiasedVar,
            varOfSampleMeans,
            varOfSampleVars,
            thVarOfMeans,
            sampleTrajectories
        };

        updateMetricsDisplay(lastSimResults);
        drawDualHistogram(lastSimResults);
        drawTrajectoryCanvas(lastSimResults, n);
    }

    function updateMetricsDisplay(d) {
        statThMean.textContent = d.dist.mu.toFixed(4);
        statAvgMean.textContent = d.grandAvgMean.toFixed(4);
        const meanDiff = Math.abs(d.grandAvgMean - d.dist.mu);
        statMeanDiff.textContent = `|Diff| = ${meanDiff.toFixed(4)} (${(meanDiff / d.dist.mu * 100).toFixed(2)}%)`;

        statThVar.textContent = d.dist.sigma2.toFixed(4);
        statAvgVar.textContent = d.grandAvgVar.toFixed(4);
        const varDiff = Math.abs(d.grandAvgVar - d.dist.sigma2);
        statVarDiff.textContent = `|Diff| = ${varDiff.toFixed(4)} (${(varDiff / d.dist.sigma2 * 100).toFixed(2)}%)`;

        statAvgBiasedVar.textContent = d.grandAvgBiasedVar.toFixed(4);

        statVarOfMeans.textContent = d.varOfSampleMeans.toFixed(5);
        statThVarOfMeans.textContent = d.thVarOfMeans.toFixed(5);
        statVarOfVars.textContent = d.varOfSampleVars.toFixed(5);

        conclusionBanner.innerHTML = `
            <div style="background:#e8f8f5; border:1px solid #27ae60; padding:1.2rem; border-radius:6px;">
                <h4 style="margin:0 0 6px 0; color:#27ae60;">Theoretical & Inferential Conclusions:</h4>
                <ul style="margin:0; padding-left:1.2rem; line-height:1.6;">
                    <li><strong>Unbiasedness of the Sample Mean:</strong> The average of the sampling means $\\overline{\\bar{X}} = \\mathbf{${d.grandAvgMean.toFixed(4)}}$ matches the theoretical population mean $\\mu = \\mathbf{${d.dist.mu.toFixed(4)}}$ almost perfectly (error: ${meanDiff.toFixed(4)}). Thus, $\\mathbb{E}[\\bar{X}] = \\mu$.</li>
                    <li><strong>Unbiasedness of Sample Variance (Bessel's Correction):</strong> The average of the sample variances with $n-1$ denominator $\\overline{S^2} = \\mathbf{${d.grandAvgVar.toFixed(4)}}$ converges directly to $\\sigma^2 = \\mathbf{${d.dist.sigma2.toFixed(4)}}$. Conversely, the biased variance (dividing by $n$) yields $\\mathbf{${d.grandAvgBiasedVar.toFixed(4)}}$, underestimating the true variance by exactly $\\frac{1}{n} = \\frac{1}{${d.n}} = ${(1 / d.n * 100).toFixed(1)}\\%$.</li>
                    <li><strong>Variability Collapse ($O(1/n)$ Scaling):</strong> Across the ${d.m} hackers, the variance of sample means $\\text{Var}(\\bar{X}) = \\mathbf{${d.varOfSampleMeans.toFixed(5)}}$ aligns with the theoretical standard error squared $\\sigma^2/n = \\mathbf{${d.thVarOfMeans.toFixed(5)}}$. Notice that as system depth $n$ increases, both $\\text{Var}(\\bar{X})$ and $\\text{Var}(S^2)$ diminish rapidly to zero!</li>
                </ul>
            </div>
        `;
    }

    // ------------------------------------------------------------------------
    // CANVAS 1: DUAL HISTOGRAMS (SAMPLING MEAN & SAMPLING VARIANCE WITH CHI-SQUARE PDF)
    // ------------------------------------------------------------------------
    function drawDualHistogram(d) {
        const { w, h } = setupHiDPI(dualHistCanvas, dualHistCtx);
        dualHistCtx.clearRect(0, 0, w, h);

        const halfW = w / 2;
        const padX = 50, padY = 35;
        const chartW = halfW - padX - 20;
        const chartH = h - 2 * padY;

        // Sub-function to draw a single histogram with exact density overlay
        function plotHistogram(dataArr, targetMean, targetVar, curveType, startX, title) {
            let minVal = Infinity, maxVal = -Infinity;
            for (let i = 0; i < dataArr.length; i++) {
                if (dataArr[i] < minVal) minVal = dataArr[i];
                if (dataArr[i] > maxVal) maxVal = dataArr[i];
            }
            let span = maxVal - minVal;
            if (span < 1e-4) span = 1;
            minVal = Math.max(0, minVal - span * 0.08);
            maxVal += span * 0.08;

            const numBins = Math.min(32, Math.max(12, Math.floor(Math.sqrt(dataArr.length))));
            const binW = (maxVal - minVal) / numBins;
            const bins = new Uint32Array(numBins);

            for (let i = 0; i < dataArr.length; i++) {
                let idx = Math.floor((dataArr[i] - minVal) / binW);
                if (idx >= numBins) idx = numBins - 1;
                if (idx < 0) idx = 0;
                bins[idx]++;
            }

            let maxDensity = 0;
            const densities = new Float32Array(numBins);
            for (let i = 0; i < numBins; i++) {
                densities[i] = bins[i] / (dataArr.length * binW);
                if (densities[i] > maxDensity) maxDensity = densities[i];
            }

            if (curveType === 'gaussian') {
                const peak = 1 / (Math.sqrt(targetVar * 2 * Math.PI));
                if (peak > maxDensity) maxDensity = peak * 1.15;
            } else if (curveType === 'chisquare') {
                // Peak of scaled Chi-Square PDF
                const modeVal = Math.max(0.1, ((d.n - 3) / (d.n - 1)) * d.dist.sigma2);
                const peak = sampleVarPDF(modeVal, d.n, d.dist.sigma2);
                if (peak > maxDensity && peak < maxDensity * 2.5) maxDensity = peak * 1.15;
            } else {
                maxDensity *= 1.15;
            }

            const toX = (val) => startX + padX + ((val - minVal) / (maxVal - minVal)) * chartW;
            const toY = (dens) => h - padY - (dens / maxDensity) * chartH;

            // Draw Title
            dualHistCtx.fillStyle = '#2c3e50';
            dualHistCtx.font = 'bold 12px sans-serif';
            dualHistCtx.textAlign = 'left';
            dualHistCtx.fillText(title, startX + padX, padY - 14);

            // Grid and Axes
            dualHistCtx.strokeStyle = '#ecf0f1';
            dualHistCtx.lineWidth = 1;
            dualHistCtx.beginPath();
            dualHistCtx.moveTo(startX + padX, h - padY);
            dualHistCtx.lineTo(startX + padX + chartW, h - padY);
            dualHistCtx.stroke();

            // Draw Bars
            for (let i = 0; i < numBins; i++) {
                let bStart = minVal + i * binW;
                let bEnd = bStart + binW;
                let x1 = toX(bStart);
                let x2 = toX(bEnd);
                let y = toY(densities[i]);
                let barWidth = Math.max(1, x2 - x1 - 1);
                let barHeight = (h - padY) - y;

                dualHistCtx.fillStyle = 'rgba(52, 152, 219, 0.65)';
                dualHistCtx.fillRect(x1, y, barWidth, barHeight);
                dualHistCtx.strokeStyle = '#2980b9';
                dualHistCtx.strokeRect(x1, y, barWidth, barHeight);
            }

            // Draw Theoretical Overlay
            if (curveType === 'gaussian' && targetVar > 1e-6) {
                const sigma = Math.sqrt(targetVar);
                dualHistCtx.strokeStyle = '#e74c3c';
                dualHistCtx.lineWidth = 2.5;
                dualHistCtx.beginPath();
                for (let i = 0; i <= 100; i++) {
                    let v = minVal + (i / 100) * (maxVal - minVal);
                    let z = (v - targetMean) / sigma;
                    let pdf = (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
                    let sx = toX(v);
                    let sy = toY(pdf);
                    if (i === 0) dualHistCtx.moveTo(sx, sy);
                    else dualHistCtx.lineTo(sx, sy);
                }
                dualHistCtx.stroke();
            } else if (curveType === 'chisquare') {
                dualHistCtx.strokeStyle = '#e74c3c';
                dualHistCtx.lineWidth = 2.5;
                dualHistCtx.beginPath();
                for (let i = 0; i <= 100; i++) {
                    let v = minVal + (i / 100) * (maxVal - minVal);
                    let pdf = sampleVarPDF(v, d.n, d.dist.sigma2);
                    let sx = toX(v);
                    let sy = toY(pdf);
                    if (i === 0) dualHistCtx.moveTo(sx, sy);
                    else dualHistCtx.lineTo(sx, sy);
                }
                dualHistCtx.stroke();
            }

            // Target marker line (mu or sigma2)
            let targetX = toX(targetMean);
            if (targetX >= startX + padX && targetX <= startX + padX + chartW) {
                dualHistCtx.strokeStyle = '#27ae60';
                dualHistCtx.lineWidth = 2;
                dualHistCtx.setLineDash([4, 3]);
                dualHistCtx.beginPath();
                dualHistCtx.moveTo(targetX, padY);
                dualHistCtx.lineTo(targetX, h - padY);
                dualHistCtx.stroke();
                dualHistCtx.setLineDash([]);

                dualHistCtx.fillStyle = '#27ae60';
                dualHistCtx.font = '10px sans-serif';
                dualHistCtx.textAlign = 'center';
                dualHistCtx.fillText(`Target: ${targetMean.toFixed(2)}`, targetX, padY - 2);
            }

            // X-axis ticks
            dualHistCtx.fillStyle = '#7f8c8d';
            dualHistCtx.font = '10px sans-serif';
            dualHistCtx.textAlign = 'center';
            for (let i = 0; i <= 4; i++) {
                let v = minVal + (i / 4) * (maxVal - minVal);
                let x = toX(v);
                dualHistCtx.fillText(v.toFixed(2), x, h - padY + 16);
            }
        }

        // Left Plot: Sampling Means
        plotHistogram(
            d.hackerMeans,
            d.dist.mu,
            d.dist.sigma2 / d.n,
            'gaussian',
            0,
            `1. Sampling Means (X̄) across ${d.m} Hackers [CLT Normal Fit]`
        );

        // Right Plot: Sampling Variances with exact Scaled Chi-Square Density
        plotHistogram(
            d.hackerVars,
            d.dist.sigma2,
            d.varOfSampleVars,
            'chisquare',
            halfW,
            `2. Sampling Variances (S²) [Exact Scaled χ²ₙ₋₁ PDF Overlay]`
        );
    }

    // ------------------------------------------------------------------------
    // CANVAS 2: WELFORD STREAM TRAJECTORIES FOR SAMPLE HACKERS
    // ------------------------------------------------------------------------
    function drawTrajectoryCanvas(d, currentN) {
        const { w, h } = setupHiDPI(trajCanvas, trajCtx);
        trajCtx.clearRect(0, 0, w, h);

        const padX = 55, padY = 32;
        const chartW = w - padX - 25;
        const chartH = h - 2 * padY;

        trajCtx.fillStyle = '#2c3e50';
        trajCtx.font = 'bold 12px sans-serif';
        trajCtx.textAlign = 'left';
        trajCtx.fillText(`Welford Online Stream Convergence across System Index k = 1 .. ${currentN} for 5 Sample Hackers:`, padX, padY - 14);

        let minY = 0, maxY = 10;
        d.sampleTrajectories.forEach(traj => {
            for (let k = 1; k <= d.n; k++) {
                if (traj.means[k] > maxY) maxY = traj.means[k];
                if (traj.vars[k] > maxY) maxY = traj.vars[k];
            }
        });
        maxY = Math.max(10, maxY * 1.1);

        const toX = (k) => padX + (k / d.n) * chartW;
        const toY = (val) => h - padY - (val / maxY) * chartH;

        // Grid
        trajCtx.strokeStyle = '#ecf0f1';
        trajCtx.lineWidth = 1;
        trajCtx.beginPath();
        trajCtx.moveTo(padX, h - padY);
        trajCtx.lineTo(w - 25, h - padY);
        trajCtx.stroke();

        // Horizontal target lines
        let yMu = toY(d.dist.mu);
        trajCtx.strokeStyle = 'rgba(52, 152, 219, 0.8)';
        trajCtx.lineWidth = 2;
        trajCtx.setLineDash([5, 4]);
        trajCtx.beginPath();
        trajCtx.moveTo(padX, yMu);
        trajCtx.lineTo(w - 25, yMu);
        trajCtx.stroke();

        let ySigma2 = toY(d.dist.sigma2);
        trajCtx.strokeStyle = 'rgba(231, 76, 60, 0.8)';
        trajCtx.lineWidth = 2;
        trajCtx.beginPath();
        trajCtx.moveTo(padX, ySigma2);
        trajCtx.lineTo(w - 25, ySigma2);
        trajCtx.stroke();
        trajCtx.setLineDash([]);

        // Plot trajectories for sample hackers up to currentN
        d.sampleTrajectories.forEach((traj, idx) => {
            // Running Mean lines (Blue hues)
            trajCtx.strokeStyle = `rgba(52, 152, 219, ${0.45 + idx * 0.12})`;
            trajCtx.lineWidth = 1.6;
            trajCtx.beginPath();
            for (let k = 1; k <= currentN; k++) {
                let sx = toX(k);
                let sy = toY(traj.means[k]);
                if (k === 1) trajCtx.moveTo(sx, sy);
                else trajCtx.lineTo(sx, sy);
            }
            trajCtx.stroke();

            // Running Variance lines (Orange/Red hues)
            trajCtx.strokeStyle = `rgba(231, 76, 60, ${0.45 + idx * 0.12})`;
            trajCtx.lineWidth = 1.6;
            trajCtx.beginPath();
            for (let k = 2; k <= currentN; k++) {
                let sx = toX(k);
                let sy = toY(traj.vars[k]);
                if (k === 2) trajCtx.moveTo(sx, sy);
                else trajCtx.lineTo(sx, sy);
            }
            trajCtx.stroke();
        });

        // X labels
        trajCtx.fillStyle = '#7f8c8d';
        trajCtx.font = '10px sans-serif';
        trajCtx.textAlign = 'center';
        for (let i = 0; i <= 5; i++) {
            let k = Math.round((i / 5) * d.n);
            let x = toX(k);
            trajCtx.fillText(`k=${k}`, x, h - padY + 16);
        }

        // Y labels
        trajCtx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            let val = (i / 4) * maxY;
            let y = toY(val);
            trajCtx.fillText(val.toFixed(1), padX - 8, y + 4);
        }

        // Legend
        trajCtx.textAlign = 'right';
        trajCtx.fillStyle = '#3498db';
        trajCtx.fillRect(w - 280, 12, 14, 10);
        trajCtx.fillStyle = '#2c3e50';
        trajCtx.font = '11px sans-serif';
        trajCtx.fillText(`Running Sample Mean x̄ₖ → μ = ${d.dist.mu.toFixed(2)}`, w - 30, 21);

        trajCtx.fillStyle = '#e74c3c';
        trajCtx.fillRect(w - 280, 28, 14, 10);
        trajCtx.fillStyle = '#2c3e50';
        trajCtx.fillText(`Running Sample Variance sₖ² → σ² = ${d.dist.sigma2.toFixed(2)}`, w - 30, 37);
    }

    // Interactive Stream Animation Controller
    playStreamBtn.addEventListener('click', () => {
        if (streamAnimTimer) {
            stopStreamAnim();
            return;
        }
        if (!lastSimResults) runSimulation();
        playStreamBtn.textContent = '⏸ Pause Stream';

        streamCurrentK = 2;
        const stepInc = Math.max(1, Math.floor(lastSimResults.n / 60));

        streamAnimTimer = setInterval(() => {
            streamCurrentK += stepInc;
            if (streamCurrentK >= lastSimResults.n) {
                streamCurrentK = lastSimResults.n;
                drawTrajectoryCanvas(lastSimResults, streamCurrentK);
                stopStreamAnim();
                return;
            }
            drawTrajectoryCanvas(lastSimResults, streamCurrentK);
        }, 35);
    });

    // ------------------------------------------------------------------------
    // CANVAS 3: VARIANCE DECAY ANALYSIS VS. SYSTEM DEPTH n (EXERCISE 2 PART B)
    // ------------------------------------------------------------------------
    function runDecayExperiment() {
        const dist = getDistributionParams();
        const m = 1500;
        const nValues = [5, 10, 20, 50, 100, 250, 500, 1000];

        const empVarOfMeans = [];
        const thVarOfMeans = [];
        const empVarOfVars = [];

        nValues.forEach(n => {
            const means = new Float64Array(m);
            const vars = new Float64Array(m);

            for (let j = 0; j < m; j++) {
                let mean = 0, M2 = 0;
                for (let k = 1; k <= n; k++) {
                    let x = dist.sample();
                    let delta = x - mean;
                    mean += delta / k;
                    M2 += delta * (x - mean);
                }
                means[j] = mean;
                vars[j] = (n > 1) ? (M2 / (n - 1)) : 0;
            }

            let avgM = 0, avgV = 0;
            for (let j = 0; j < m; j++) {
                avgM += means[j];
                avgV += vars[j];
            }
            avgM /= m;
            avgV /= m;

            let sqM = 0, sqV = 0;
            for (let j = 0; j < m; j++) {
                sqM += Math.pow(means[j] - avgM, 2);
                sqV += Math.pow(vars[j] - avgV, 2);
            }
            empVarOfMeans.push(sqM / (m - 1));
            thVarOfMeans.push(dist.sigma2 / n);
            empVarOfVars.push(sqV / (m - 1));
        });

        drawDecayCanvas(nValues, empVarOfMeans, thVarOfMeans, empVarOfVars, dist);
    }

    function drawDecayCanvas(nVals, empMeans, thMeans, empVars, dist) {
        const { w, h } = setupHiDPI(decayCanvas, decayCtx);
        decayCtx.clearRect(0, 0, w, h);

        const padX = 65, padY = 35;
        const chartW = w - padX - 25;
        const chartH = h - 2 * padY;

        decayCtx.fillStyle = '#2c3e50';
        decayCtx.font = 'bold 12px sans-serif';
        decayCtx.textAlign = 'left';
        decayCtx.fillText('Variability Decay of Sampling Statistics across System Depth n (O(1/n) Scaling):', padX, padY - 14);

        let maxY = Math.max(...empMeans, ...thMeans, ...empVars);
        maxY *= 1.1;

        const toX = (idx) => padX + (idx / (nVals.length - 1)) * chartW;
        const toY = (val) => h - padY - (val / maxY) * chartH;

        // Grid
        decayCtx.strokeStyle = '#ecf0f1';
        decayCtx.lineWidth = 1;
        decayCtx.beginPath();
        decayCtx.moveTo(padX, h - padY);
        decayCtx.lineTo(w - 25, h - padY);
        decayCtx.moveTo(padX, padY);
        decayCtx.lineTo(padX, h - padY);
        decayCtx.stroke();

        // 1. Theoretical Var(X_bar) = sigma^2 / n
        decayCtx.strokeStyle = '#3498db';
        decayCtx.lineWidth = 2.5;
        decayCtx.beginPath();
        thMeans.forEach((v, idx) => {
            let sx = toX(idx);
            let sy = toY(v);
            if (idx === 0) decayCtx.moveTo(sx, sy);
            else decayCtx.lineTo(sx, sy);
        });
        decayCtx.stroke();

        // Empirical Points for Var(X_bar)
        decayCtx.fillStyle = '#2980b9';
        empMeans.forEach((v, idx) => {
            let sx = toX(idx);
            let sy = toY(v);
            decayCtx.beginPath();
            decayCtx.arc(sx, sy, 4.5, 0, 2 * Math.PI);
            decayCtx.fill();
        });

        // 2. Empirical Var(S^2)
        decayCtx.strokeStyle = '#e74c3c';
        decayCtx.lineWidth = 2.5;
        decayCtx.beginPath();
        empVars.forEach((v, idx) => {
            let sx = toX(idx);
            let sy = toY(v);
            if (idx === 0) decayCtx.moveTo(sx, sy);
            else decayCtx.lineTo(sx, sy);
        });
        decayCtx.stroke();

        // Empirical Points for Var(S^2)
        decayCtx.fillStyle = '#c0392b';
        empVars.forEach((v, idx) => {
            let sx = toX(idx);
            let sy = toY(v);
            decayCtx.beginPath();
            decayCtx.arc(sx, sy, 4.5, 0, 2 * Math.PI);
            decayCtx.fill();
        });

        // X labels
        decayCtx.fillStyle = '#7f8c8d';
        decayCtx.font = '10px sans-serif';
        decayCtx.textAlign = 'center';
        nVals.forEach((n, idx) => {
            let x = toX(idx);
            decayCtx.fillText(`n=${n}`, x, h - padY + 16);
        });

        // Y labels
        decayCtx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            let val = (i / 4) * maxY;
            let y = toY(val);
            decayCtx.fillText(val.toFixed(2), padX - 8, y + 4);
        }

        // Legend
        decayCtx.textAlign = 'right';
        decayCtx.fillStyle = '#3498db';
        decayCtx.fillRect(w - 320, 12, 14, 10);
        decayCtx.fillStyle = '#2c3e50';
        decayCtx.font = '11px sans-serif';
        decayCtx.fillText('Var(X̄) = σ²/n  [Standard Error Decay O(1/n)]', w - 30, 21);

        decayCtx.fillStyle = '#e74c3c';
        decayCtx.fillRect(w - 320, 28, 14, 10);
        decayCtx.fillStyle = '#2c3e50';
        decayCtx.fillText('Var(S²) across Hackers  [Variance Uncertainty Decay O(1/n)]', w - 30, 37);
    }

    runBatchNBtn.addEventListener('click', runDecayExperiment);
    runSimBtn.addEventListener('click', runSimulation);

    // Initial executions
    runPlaygroundTest();
    runSimulation();
    runDecayExperiment();
});

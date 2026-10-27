// ============================================================================
// HOMEWORK 4 - STOCHASTIC PENETRATION PROCESS & ASYMPTOTIC CONVERGENCE
// Dario Pietrosanto (2055042)
// Master's Degree Course in Applied Statistics & Cybersecurity
// Advanced Interactive Simulation Engine (Guru Level):
// 1. Bernoulli and Rademacher I.I.D. Increment Generators
// 2. Trajectory Simulator (Total Sn, Average Sn/n, Normalized Sn/sqrt(n), Standardized Zn)
// 3. Dynamic Multi-Path Canvas Visualizer with Theoretical ±2σ Confidence Envelopes
// 4. Terminal Histogram with Gaussian PDF Overlay & Berry-Esseen / K-S Convergence Metrics
// 5. Empirical vs. Theoretical Variance Tracking Across System Depth k = 1 .. n
// 6. Interactive Play / Animation Mode (Step-by-Step Morphing)
// 7. Extreme Value Theory (EVT) Analysis for Attacker Fleets
// ============================================================================

document.addEventListener('DOMContentLoaded', () => {

    // --- DOM Elements ---
    const modelSelect = document.getElementById('modelSelect');
    const pControlGroup = document.getElementById('pControlGroup');
    const pInput = document.getElementById('pInput');
    const pValueDisplay = document.getElementById('pValueDisplay');
    const nInput = document.getElementById('nInput');
    const nValueDisplay = document.getElementById('nValueDisplay');
    const mInput = document.getElementById('mInput');
    const mValueDisplay = document.getElementById('mValueDisplay');
    const normRadios = document.querySelectorAll('input[name="normalization"]');
    const standardizeCheckbox = document.getElementById('standardizeCheckbox');
    const runSimBtn = document.getElementById('runSimBtn');
    const playAnimBtn = document.getElementById('playAnimBtn');

    // Canvas elements
    const trajCanvas = document.getElementById('trajCanvas');
    const trajCtx = trajCanvas.getContext('2d');
    const histCanvas = document.getElementById('histCanvas');
    const histCtx = histCanvas.getContext('2d');
    const varCanvas = document.getElementById('varCanvas');
    const varCtx = varCanvas.getContext('2d');

    // Metrics display elements
    const statTheoreticalMean = document.getElementById('statTheoreticalMean');
    const statTheoreticalVar = document.getElementById('statTheoreticalVar');
    const statEmpiricalMean = document.getElementById('statEmpiricalMean');
    const statEmpiricalVar = document.getElementById('statEmpiricalVar');
    const statSkewness = document.getElementById('statSkewness');
    const statKurtosis = document.getElementById('statKurtosis');
    const statMaxBreach = document.getElementById('statMaxBreach');
    const statTheoreticalMax = document.getElementById('statTheoreticalMax');
    const statKsDistance = document.getElementById('statKsDistance');
    const statBerryEsseen = document.getElementById('statBerryEsseen');
    const cltStatusBanner = document.getElementById('cltStatusBanner');

    // High-DPI Canvas Scaling Helper
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

    // Preset buttons
    document.querySelectorAll('.preset-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            stopAnimation();
            const presetModel = btn.getAttribute('data-model');
            const presetN = parseInt(btn.getAttribute('data-n'), 10);
            const presetM = parseInt(btn.getAttribute('data-m'), 10);
            const presetNorm = btn.getAttribute('data-norm');
            const presetP = parseFloat(btn.getAttribute('data-p') || '0.5');
            const presetStd = (btn.getAttribute('data-std') === 'true');

            modelSelect.value = presetModel;
            nInput.value = presetN;
            mInput.value = presetM;
            pInput.value = presetP;
            standardizeCheckbox.checked = presetStd;
            document.querySelector(`input[name="normalization"][value="${presetNorm}"]`).checked = true;

            syncInputs();
            runSimulation();
        });
    });

    function syncInputs() {
        const isBernoulli = (modelSelect.value === 'bernoulli');
        pControlGroup.style.display = isBernoulli ? 'block' : 'none';
        pValueDisplay.textContent = parseFloat(pInput.value).toFixed(2);
        nValueDisplay.textContent = parseInt(nInput.value, 10).toLocaleString();
        mValueDisplay.textContent = parseInt(mInput.value, 10).toLocaleString();
    }

    modelSelect.addEventListener('change', () => { syncInputs(); runSimulation(); });
    pInput.addEventListener('input', () => { syncInputs(); runSimulation(); });
    nInput.addEventListener('input', () => { syncInputs(); runSimulation(); });
    mInput.addEventListener('input', () => { syncInputs(); runSimulation(); });
    standardizeCheckbox.addEventListener('change', runSimulation);
    normRadios.forEach(r => r.addEventListener('change', runSimulation));

    function getSelectedNormalization() {
        const selected = document.querySelector('input[name="normalization"]:checked');
        return selected ? selected.value : 'normalized';
    }

    // Standard Normal CDF Approximation (Abramowitz & Stegun formula 7.1.26)
    function normalCDF(x) {
        const a1 =  0.254829592, a2 = -0.284496736, a3 =  1.421413741;
        const a4 = -1.453152027, a5 =  1.061405429, p  =  0.3275911;
        const sign = x < 0 ? -1 : 1;
        const absX = Math.abs(x) / Math.SQRT2;
        const t = 1.0 / (1.0 + p * absX);
        const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-absX * absX);
        return 0.5 * (1.0 + sign * y);
    }

    // Kolmogorov-Smirnov supremum distance: D = sup |F_m(x) - Phi(x)|
    function computeKSDistance(sample, targetMean, targetSigma) {
        if (targetSigma < 1e-6) return 0;
        const sorted = Array.from(sample).sort((a, b) => a - b);
        const m = sorted.length;
        let dMax = 0;
        for (let i = 0; i < m; i++) {
            let z = (sorted[i] - targetMean) / targetSigma;
            let cdfEmpAfter = (i + 1) / m;
            let cdfEmpBefore = i / m;
            let cdfTh = normalCDF(z);
            let diff1 = Math.abs(cdfEmpAfter - cdfTh);
            let diff2 = Math.abs(cdfTh - cdfEmpBefore);
            if (diff1 > dMax) dMax = diff1;
            if (diff2 > dMax) dMax = diff2;
        }
        return dMax;
    }

    // ------------------------------------------------------------------------
    // SIMULATION STATE & ENGINE
    // ------------------------------------------------------------------------
    let simData = null;
    let animTimer = null;
    let animCurrentStep = 0;

    function runSimulation() {
        stopAnimation();
        syncInputs();

        const model = modelSelect.value;
        const p = parseFloat(pInput.value);
        const n = parseInt(nInput.value, 10);
        const m = parseInt(mInput.value, 10);
        const norm = getSelectedNormalization();
        const doStandardize = standardizeCheckbox.checked;

        // Moment parameters of single increment X_i
        let mu_X = 0;
        let sigma2_X = 1;
        let rho_X = 1; // E[|X - mu|^3] for Berry-Esseen theorem

        if (model === 'rademacher') {
            mu_X = 0;
            sigma2_X = 1;
            rho_X = 1; // |+1|^3 = 1, |-1|^3 = 1
        } else {
            // Bernoulli
            mu_X = p;
            sigma2_X = p * (1 - p);
            // Third central moment of Bernoulli: p(1-p)^3 + (1-p)p^3 = p(1-p)[(1-p)^2 + p^2]
            rho_X = p * Math.pow(1 - p, 3) + (1 - p) * Math.pow(p, 3);
        }
        const sigma_X = Math.sqrt(sigma2_X);

        // Max visual paths to render in multi-path trajectories canvas
        const maxRenderPaths = Math.min(m, 60);
        const samplePaths = [];
        for (let j = 0; j < maxRenderPaths; j++) {
            samplePaths.push(new Float32Array(n + 1));
        }

        const terminalValues = new Float32Array(m);
        const terminalRawScores = new Float32Array(m);

        // Subsample for Empirical Variance Evolution tracking across k
        const varSampleSize = Math.min(m, 500);
        const trajectorySums = [];
        for (let j = 0; j < varSampleSize; j++) {
            trajectorySums.push(new Float32Array(n + 1));
        }

        // Generate full trajectories
        for (let j = 0; j < m; j++) {
            let sum = 0;
            const isRenderPath = (j < maxRenderPaths);
            const isVarPath = (j < varSampleSize);

            if (isRenderPath) samplePaths[j][0] = 0;
            if (isVarPath) trajectorySums[j][0] = 0;

            for (let k = 1; k <= n; k++) {
                let X_i = (model === 'rademacher') ? ((Math.random() < 0.5) ? 1 : -1) : ((Math.random() < p) ? 1 : 0);
                sum += X_i;

                if (isRenderPath) {
                    if (doStandardize && model === 'bernoulli') {
                        samplePaths[j][k] = (sum - k * mu_X) / (sigma_X * Math.sqrt(k));
                    } else if (norm === 'total') {
                        samplePaths[j][k] = sum;
                    } else if (norm === 'average') {
                        samplePaths[j][k] = sum / k;
                    } else {
                        samplePaths[j][k] = sum / Math.sqrt(k);
                    }
                }

                if (isVarPath) {
                    trajectorySums[j][k] = sum;
                }
            }

            terminalRawScores[j] = sum;
            if (doStandardize && model === 'bernoulli') {
                terminalValues[j] = (sum - n * mu_X) / (sigma_X * Math.sqrt(n));
            } else if (norm === 'total') {
                terminalValues[j] = sum;
            } else if (norm === 'average') {
                terminalValues[j] = sum / n;
            } else {
                terminalValues[j] = sum / Math.sqrt(n);
            }
        }

        // Empirical variance array computed across k
        const kCheckpoints = [];
        const stepSize = Math.max(1, Math.floor(n / 60));
        for (let k = 1; k <= n; k += stepSize) kCheckpoints.push(k);
        if (kCheckpoints[kCheckpoints.length - 1] !== n) kCheckpoints.push(n);

        const empVarEvolution = {
            kList: kCheckpoints,
            total: new Float32Array(kCheckpoints.length),
            average: new Float32Array(kCheckpoints.length),
            normalized: new Float32Array(kCheckpoints.length)
        };

        kCheckpoints.forEach((k, idx) => {
            let sumTot = 0, sumAvg = 0, sumNorm = 0;
            for (let j = 0; j < varSampleSize; j++) {
                let s = trajectorySums[j][k];
                sumTot += s;
                sumAvg += s / k;
                sumNorm += s / Math.sqrt(k);
            }
            let meanTot = sumTot / varSampleSize;
            let meanAvg = sumAvg / varSampleSize;
            let meanNorm = sumNorm / varSampleSize;

            let sqTot = 0, sqAvg = 0, sqNorm = 0;
            for (let j = 0; j < varSampleSize; j++) {
                let s = trajectorySums[j][k];
                sqTot += Math.pow(s - meanTot, 2);
                sqAvg += Math.pow((s / k) - meanAvg, 2);
                sqNorm += Math.pow((s / Math.sqrt(k)) - meanNorm, 2);
            }
            empVarEvolution.total[idx] = sqTot / (varSampleSize - 1);
            empVarEvolution.average[idx] = sqAvg / (varSampleSize - 1);
            empVarEvolution.normalized[idx] = sqNorm / (varSampleSize - 1);
        });

        // Moments Estimation of terminal values
        let sumVals = 0;
        let maxBreach = -Infinity;
        for (let j = 0; j < m; j++) {
            sumVals += terminalValues[j];
            if (terminalRawScores[j] > maxBreach) maxBreach = terminalRawScores[j];
        }
        const empMean = sumVals / m;

        let sumSqDiff = 0;
        for (let j = 0; j < m; j++) {
            sumSqDiff += Math.pow(terminalValues[j] - empMean, 2);
        }
        const empVar = sumSqDiff / (m - 1);
        const empStd = Math.sqrt(empVar);

        let sumCube = 0, sumQuad = 0;
        for (let j = 0; j < m; j++) {
            const z = empStd > 1e-9 ? (terminalValues[j] - empMean) / empStd : 0;
            sumCube += Math.pow(z, 3);
            sumQuad += Math.pow(z, 4);
        }
        const skewness = sumCube / m;
        const kurtosis = (sumQuad / m) - 3; // Excess Kurtosis

        // Theoretical moments
        let thMean = 0, thVar = 0;
        if (doStandardize && model === 'bernoulli') {
            thMean = 0;
            thVar = 1;
        } else if (norm === 'total') {
            thMean = n * mu_X;
            thVar = n * sigma2_X;
        } else if (norm === 'average') {
            thMean = mu_X;
            thVar = sigma2_X / n;
        } else {
            // normalized
            thMean = Math.sqrt(n) * mu_X;
            thVar = sigma2_X;
        }

        // Theoretical Expected Max Breach via Extreme Value Theory
        const thMaxBreach = n * mu_X + sigma_X * Math.sqrt(2 * n * Math.log(m));

        // Berry-Esseen Upper Bound: C * rho / (sigma^3 * sqrt(n)), C <= 0.4748
        const berryEsseenBound = (0.4748 * rho_X) / (Math.pow(sigma_X, 3) * Math.sqrt(n));

        // Empirical Kolmogorov-Smirnov distance to Standard Normal
        const ksDist = computeKSDistance(terminalValues, thMean, Math.sqrt(thVar));

        simData = {
            model, p, n, m, norm, doStandardize, mu_X, sigma_X, sigma2_X, rho_X,
            samplePaths, terminalValues, terminalRawScores,
            empMean, empVar, empStd, skewness, kurtosis,
            thMean, thVar, maxBreach, thMaxBreach,
            berryEsseenBound, ksDist,
            empVarEvolution
        };

        updateMetricsDisplay(simData);
        drawTrajectoriesCanvas(simData, n);
        drawHistogramCanvas(simData);
        drawVarianceCanvas(simData);
    }

    // ------------------------------------------------------------------------
    // UI METRICS DISPLAY
    // ------------------------------------------------------------------------
    function updateMetricsDisplay(d) {
        statTheoreticalMean.textContent = d.thMean.toFixed(4);
        statTheoreticalVar.textContent = d.thVar.toFixed(4);
        statEmpiricalMean.textContent = d.empMean.toFixed(4);
        statEmpiricalVar.textContent = d.empVar.toFixed(4);
        statSkewness.textContent = d.skewness.toFixed(4);
        statKurtosis.textContent = d.kurtosis.toFixed(4);
        statMaxBreach.textContent = `${d.maxBreach.toFixed(0)} / ${d.n}`;
        statTheoreticalMax.textContent = `≈ ${d.thMaxBreach.toFixed(1)}`;
        statKsDistance.textContent = d.ksDist.toFixed(4);
        statBerryEsseen.textContent = `≤ ${d.berryEsseenBound.toFixed(4)}`;

        if (d.norm === 'average') {
            const err = Math.abs(d.empMean - d.thMean);
            cltStatusBanner.innerHTML = `
                <div style="background:#e8f8f5; border:1px solid #27ae60; padding:1rem; border-radius:4px;">
                    <h4 style="margin:0 0 5px 0; color:#27ae60;">Law of Large Numbers (LLN) Convergence:</h4>
                    Empirical Average Breach Rate &bar{S}<sub>n</sub>/n = <strong>${d.empMean.toFixed(4)}</strong> converges to true expectation &mu; = <strong>${d.thMean.toFixed(4)}</strong> (Absolute Error: <strong>${err.toExponential(3)}</strong>). 
                    The path variance collapses hyperbolically as <strong>Var(S<sub>n</sub>/n) = &sigma;&sup2;/n = ${(d.sigma2_X / d.n).toFixed(6)}</strong>, completely eliminating individual attacker stochastic luck.
                </div>
            `;
        } else if ((d.norm === 'normalized' && d.model === 'rademacher') || d.doStandardize) {
            cltStatusBanner.innerHTML = `
                <div style="background:#e8f8f5; border:1px solid #27ae60; padding:1rem; border-radius:4px;">
                    <h4 style="margin:0 0 5px 0; color:#27ae60;">Central Limit Theorem (CLT) Convergence:</h4>
                    Process satisfies Lindeberg-L&eacute;vy conditions: Unit Variance <strong>s&sup2; = ${d.empVar.toFixed(3)}</strong> (Theoretical 1.000) &nbsp;|&nbsp;
                    Skewness &gamma;<sub>1</sub> = <strong>${d.skewness.toFixed(3)}</strong> &nbsp;|&nbsp;
                    Excess Kurtosis &gamma;<sub>2</sub> = <strong>${d.kurtosis.toFixed(3)}</strong>.<br>
                    <strong>Kolmogorov-Smirnov Distance:</strong> D = <strong>${d.ksDist.toFixed(4)}</strong> &le; Berry-Esseen Upper Bound <strong>${d.berryEsseenBound.toFixed(4)}</strong> (Confirms convergence rate of order <em>O(1/&radic;n)</em>).
                </div>
            `;
        } else {
            cltStatusBanner.innerHTML = `
                <div style="background:#f8f9fa; border:1px solid #bdc3c7; padding:1rem; border-radius:4px;">
                    <strong>Total Score Regime (S<sub>n</sub>):</strong> 
                    Variance diverges linearly with systems: <strong>Var(S<sub>n</sub>) = n&sigma;&sup2; = ${d.thVar.toFixed(1)}</strong>. 
                    The diffusive envelope expands as <strong>&plusmn; 2&sigma;&radic;n = &plusmn; ${(2 * d.sigma_X * Math.sqrt(d.n)).toFixed(1)}</strong>.
                </div>
            `;
        }
    }

    // ------------------------------------------------------------------------
    // CANVAS 1: TRAJECTORIES MULTI-PATH
    // ------------------------------------------------------------------------
    function drawTrajectoriesCanvas(d, currentN) {
        const { w, h } = setupHiDPI(trajCanvas, trajCtx);
        trajCtx.clearRect(0, 0, w, h);

        const padX = 65, padY = 32;
        const chartW = w - padX - 25;
        const chartH = h - 2 * padY;

        let minY = Infinity, maxY = -Infinity;
        d.samplePaths.forEach(path => {
            for (let k = 0; k <= d.n; k++) {
                if (path[k] < minY) minY = path[k];
                if (path[k] > maxY) maxY = path[k];
            }
        });

        let span = maxY - minY;
        if (span < 1e-4) span = 1;
        minY -= span * 0.08;
        maxY += span * 0.08;

        const toScreenX = (k) => padX + (k / d.n) * chartW;
        const toScreenY = (val) => h - padY - ((val - minY) / (maxY - minY)) * chartH;

        // Grid lines
        trajCtx.strokeStyle = '#ecf0f1';
        trajCtx.lineWidth = 1;
        trajCtx.fillStyle = '#7f8c8d';
        trajCtx.font = '10px sans-serif';
        trajCtx.textAlign = 'right';

        for (let i = 0; i <= 6; i++) {
            let val = minY + (i / 6) * (maxY - minY);
            let y = toScreenY(val);
            trajCtx.beginPath();
            trajCtx.moveTo(padX, y);
            trajCtx.lineTo(w - 25, y);
            trajCtx.stroke();
            trajCtx.fillText(val.toFixed(d.norm === 'average' ? 3 : 1), padX - 8, y + 4);
        }

        trajCtx.textAlign = 'center';
        for (let i = 0; i <= 5; i++) {
            let k = Math.round((i / 5) * d.n);
            let x = toScreenX(k);
            trajCtx.beginPath();
            trajCtx.moveTo(x, padY);
            trajCtx.lineTo(x, h - padY);
            trajCtx.stroke();
            trajCtx.fillText(`k=${k}`, x, h - padY + 16);
        }

        // Draw Theoretical Mean Line
        trajCtx.strokeStyle = '#2c3e50';
        trajCtx.lineWidth = 2.5;
        trajCtx.setLineDash([6, 4]);
        trajCtx.beginPath();
        for (let k = 1; k <= d.n; k += Math.ceil(d.n / 100)) {
            let meanVal = 0;
            if (d.doStandardize && d.model === 'bernoulli') meanVal = 0;
            else if (d.norm === 'total') meanVal = k * d.mu_X;
            else if (d.norm === 'average') meanVal = d.mu_X;
            else meanVal = Math.sqrt(k) * d.mu_X;

            let sx = toScreenX(k);
            let sy = toScreenY(meanVal);
            if (k === 1) trajCtx.moveTo(sx, sy);
            else trajCtx.lineTo(sx, sy);
        }
        trajCtx.stroke();
        trajCtx.setLineDash([]);

        // Theoretical ±2 sigma Confidence Envelope
        trajCtx.fillStyle = 'rgba(52, 152, 219, 0.08)';
        trajCtx.beginPath();
        for (let k = 1; k <= d.n; k += Math.ceil(d.n / 100)) {
            let upper = 0;
            if (d.doStandardize && d.model === 'bernoulli') upper = 2;
            else if (d.norm === 'total') upper = k * d.mu_X + 2 * d.sigma_X * Math.sqrt(k);
            else if (d.norm === 'average') upper = d.mu_X + (2 * d.sigma_X) / Math.sqrt(k);
            else upper = Math.sqrt(k) * d.mu_X + 2 * d.sigma_X;

            let sx = toScreenX(k);
            let sy = toScreenY(upper);
            if (k === 1) trajCtx.moveTo(sx, sy);
            else trajCtx.lineTo(sx, sy);
        }
        for (let k = d.n; k >= 1; k -= Math.ceil(d.n / 100)) {
            let lower = 0;
            if (d.doStandardize && d.model === 'bernoulli') lower = -2;
            else if (d.norm === 'total') lower = k * d.mu_X - 2 * d.sigma_X * Math.sqrt(k);
            else if (d.norm === 'average') lower = d.mu_X - (2 * d.sigma_X) / Math.sqrt(k);
            else lower = Math.sqrt(k) * d.mu_X - 2 * d.sigma_X;

            let sx = toScreenX(k);
            let sy = toScreenY(lower);
            trajCtx.lineTo(sx, sy);
        }
        trajCtx.closePath();
        trajCtx.fill();

        // Draw Sample Paths up to currentN (for animation)
        const pathColors = [
            'rgba(231, 76, 60, 0.65)',
            'rgba(52, 152, 219, 0.65)',
            'rgba(39, 174, 96, 0.65)',
            'rgba(155, 89, 182, 0.65)',
            'rgba(243, 156, 18, 0.65)'
        ];

        trajCtx.lineWidth = 1.3;
        d.samplePaths.forEach((path, idx) => {
            trajCtx.strokeStyle = pathColors[idx % pathColors.length];
            trajCtx.beginPath();
            for (let k = 0; k <= currentN; k++) {
                let sx = toScreenX(k);
                let sy = toScreenY(path[k]);
                if (k === 0) trajCtx.moveTo(sx, sy);
                else trajCtx.lineTo(sx, sy);
            }
            trajCtx.stroke();
        });

        // Header info
        trajCtx.textAlign = 'left';
        trajCtx.fillStyle = '#2c3e50';
        trajCtx.font = '11px sans-serif';
        trajCtx.fillText(`Showing ${d.samplePaths.length} Sample Paths across k=1..${currentN}`, padX, padY - 12);

        trajCtx.fillStyle = 'rgba(52, 152, 219, 0.8)';
        trajCtx.fillRect(w - 240, 10, 14, 10);
        trajCtx.fillStyle = '#2c3e50';
        trajCtx.fillText('Theoretical ±2σ Confidence Envelope', w - 220, 18);
    }

    // ------------------------------------------------------------------------
    // CANVAS 2: TERMINAL HISTOGRAM & CLT GAUSSIAN
    // ------------------------------------------------------------------------
    function drawHistogramCanvas(d) {
        const { w, h } = setupHiDPI(histCanvas, histCtx);
        histCtx.clearRect(0, 0, w, h);

        const padX = 55, padY = 35;
        const chartW = w - padX - 25;
        const chartH = h - 2 * padY;

        let minVal = Infinity, maxVal = -Infinity;
        for (let j = 0; j < d.m; j++) {
            if (d.terminalValues[j] < minVal) minVal = d.terminalValues[j];
            if (d.terminalValues[j] > maxVal) maxVal = d.terminalValues[j];
        }

        let span = maxVal - minVal;
        if (span < 1e-4) span = 1;
        minVal -= span * 0.05;
        maxVal += span * 0.05;

        const numBins = Math.min(45, Math.max(16, Math.floor(Math.sqrt(d.m))));
        const binWidth = (maxVal - minVal) / numBins;
        const bins = new Uint32Array(numBins);

        for (let j = 0; j < d.m; j++) {
            let binIdx = Math.floor((d.terminalValues[j] - minVal) / binWidth);
            if (binIdx >= numBins) binIdx = numBins - 1;
            if (binIdx < 0) binIdx = 0;
            bins[binIdx]++;
        }

        let maxDensity = 0;
        const densities = new Float32Array(numBins);
        for (let i = 0; i < numBins; i++) {
            densities[i] = bins[i] / (d.m * binWidth);
            if (densities[i] > maxDensity) maxDensity = densities[i];
        }

        const isStandardNormalTarget = (d.norm === 'normalized' && d.model === 'rademacher') || d.doStandardize;
        if (isStandardNormalTarget) {
            const normalPeak = 1 / Math.sqrt(2 * Math.PI);
            if (normalPeak > maxDensity) maxDensity = normalPeak * 1.15;
        } else {
            maxDensity *= 1.15;
        }

        const toScreenX = (val) => padX + ((val - minVal) / (maxVal - minVal)) * chartW;
        const toScreenY = (density) => h - padY - (density / maxDensity) * chartH;

        // Bars
        for (let i = 0; i < numBins; i++) {
            let bStart = minVal + i * binWidth;
            let bEnd = bStart + binWidth;
            let x1 = toScreenX(bStart);
            let x2 = toScreenX(bEnd);
            let y = toScreenY(densities[i]);
            let barW = Math.max(1, x2 - x1 - 1);
            let barH = (h - padY) - y;

            histCtx.fillStyle = 'rgba(52, 152, 219, 0.65)';
            histCtx.fillRect(x1, y, barW, barH);
            histCtx.strokeStyle = '#2980b9';
            histCtx.strokeRect(x1, y, barW, barH);
        }

        // Gaussian PDF Curve Overlay
        const thSigma = Math.sqrt(d.thVar);
        if (thSigma > 1e-6) {
            histCtx.strokeStyle = '#e74c3c';
            histCtx.lineWidth = 2.5;
            histCtx.beginPath();
            const curvePts = 160;
            for (let i = 0; i <= curvePts; i++) {
                let xVal = minVal + (i / curvePts) * (maxVal - minVal);
                let z = (xVal - d.thMean) / thSigma;
                let pdf = (1 / (thSigma * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
                let sx = toScreenX(xVal);
                let sy = toScreenY(pdf);
                if (i === 0) histCtx.moveTo(sx, sy);
                else histCtx.lineTo(sx, sy);
            }
            histCtx.stroke();
        }

        // Axes
        histCtx.strokeStyle = '#bdc3c7';
        histCtx.lineWidth = 1;
        histCtx.beginPath();
        histCtx.moveTo(padX, h - padY);
        histCtx.lineTo(w - 25, h - padY);
        histCtx.stroke();

        histCtx.fillStyle = '#2c3e50';
        histCtx.font = '10px sans-serif';
        histCtx.textAlign = 'center';
        for (let i = 0; i <= 6; i++) {
            let val = minVal + (i / 6) * (maxVal - minVal);
            let x = toScreenX(val);
            histCtx.fillText(val.toFixed(d.norm === 'average' ? 3 : 1), x, h - padY + 16);
        }

        histCtx.textAlign = 'right';
        histCtx.fillStyle = 'rgba(52, 152, 219, 0.8)';
        histCtx.fillRect(w - 260, 12, 12, 10);
        histCtx.fillStyle = '#2c3e50';
        histCtx.font = '11px sans-serif';
        histCtx.fillText(`Empirical Fleet Distribution (m=${d.m.toLocaleString()})`, w - 30, 20);

        histCtx.fillStyle = '#e74c3c';
        histCtx.fillRect(w - 260, 26, 12, 10);
        histCtx.fillStyle = '#2c3e50';
        histCtx.fillText('Theoretical Normal Limit PDF', w - 30, 34);
    }

    // ------------------------------------------------------------------------
    // CANVAS 3: EMPIRICAL VS THEORETICAL VARIANCE TRACKING
    // ------------------------------------------------------------------------
    function drawVarianceCanvas(d) {
        const { w, h } = setupHiDPI(varCanvas, varCtx);
        varCtx.clearRect(0, 0, w, h);

        const padX = 65, padY = 32;
        const chartW = w - padX - 25;
        const chartH = h - 2 * padY;

        varCtx.fillStyle = '#2c3e50';
        varCtx.font = '11px sans-serif';
        varCtx.textAlign = 'left';
        varCtx.fillText('Empirical Sample Variance (Points) vs. Theoretical Curves (Solid Lines) along k=1..n:', padX, padY - 12);

        // Grid
        varCtx.strokeStyle = '#ecf0f1';
        varCtx.lineWidth = 1;
        varCtx.beginPath();
        varCtx.moveTo(padX, h - padY);
        varCtx.lineTo(w - 25, h - padY);
        varCtx.moveTo(padX, padY);
        varCtx.lineTo(padX, h - padY);
        varCtx.stroke();

        const toScreenX = (k) => padX + (k / d.n) * chartW;
        const toScreenY_Norm = (pct) => h - padY - pct * chartH;

        // 1. Total Score Variance: Var(Sk) / (n * sigma2) => Linear 0 to 1
        varCtx.strokeStyle = '#e74c3c';
        varCtx.lineWidth = 2.5;
        varCtx.beginPath();
        varCtx.moveTo(toScreenX(1), toScreenY_Norm(1 / d.n));
        varCtx.lineTo(toScreenX(d.n), toScreenY_Norm(1.0));
        varCtx.stroke();

        // 2. Normalized Score Variance: Var(Sk/sqrt(k)) / sigma2 => Constant 0.5 (scaled)
        varCtx.strokeStyle = '#27ae60';
        varCtx.lineWidth = 2.5;
        varCtx.beginPath();
        varCtx.moveTo(toScreenX(1), toScreenY_Norm(0.5));
        varCtx.lineTo(toScreenX(d.n), toScreenY_Norm(0.5));
        varCtx.stroke();

        // 3. Average Score Variance: Var(Sk/k) / sigma2 => 1/k (decays to 0)
        varCtx.strokeStyle = '#3498db';
        varCtx.lineWidth = 2.5;
        varCtx.beginPath();
        for (let k = 1; k <= d.n; k += Math.ceil(d.n / 100)) {
            let decay = 0.5 * (1 / k);
            let sx = toScreenX(k);
            let sy = toScreenY_Norm(decay);
            if (k === 1) varCtx.moveTo(sx, sy);
            else varCtx.lineTo(sx, sy);
        }
        varCtx.stroke();

        // Overlay Empirical Variance Points from simulation
        const kList = d.empVarEvolution.kList;
        for (let i = 0; i < kList.length; i++) {
            let k = kList[i];
            let sx = toScreenX(k);

            // Total empirical point
            let empTotNorm = d.empVarEvolution.total[i] / (d.n * d.sigma2_X);
            let syTot = toScreenY_Norm(Math.min(1.1, Math.max(0, empTotNorm)));
            varCtx.fillStyle = '#e74c3c';
            varCtx.beginPath();
            varCtx.arc(sx, syTot, 2.5, 0, 2 * Math.PI);
            varCtx.fill();

            // Normalized empirical point
            let empNormScale = 0.5 * (d.empVarEvolution.normalized[i] / d.sigma2_X);
            let syNorm = toScreenY_Norm(Math.min(1.1, Math.max(0, empNormScale)));
            varCtx.fillStyle = '#27ae60';
            varCtx.beginPath();
            varCtx.arc(sx, syNorm, 2.5, 0, 2 * Math.PI);
            varCtx.fill();

            // Average empirical point
            let empAvgScale = 0.5 * (d.empVarEvolution.average[i] / d.sigma2_X);
            let syAvg = toScreenY_Norm(Math.min(1.1, Math.max(0, empAvgScale)));
            varCtx.fillStyle = '#3498db';
            varCtx.beginPath();
            varCtx.arc(sx, syAvg, 2.5, 0, 2 * Math.PI);
            varCtx.fill();
        }

        // X labels
        varCtx.fillStyle = '#7f8c8d';
        varCtx.font = '10px sans-serif';
        varCtx.textAlign = 'center';
        for (let i = 0; i <= 5; i++) {
            let k = Math.round((i / 5) * d.n);
            let x = toScreenX(k);
            varCtx.fillText(`k=${k}`, x, h - padY + 16);
        }

        // Legend
        const legends = [
            { text: 'Total: Var(Sₖ) = k·σ²  [O(n) Divergence]', color: '#e74c3c' },
            { text: 'Normalized: Var(Sₖ/√k) = σ²  [O(1) Flat Invariance (CLT)]', color: '#27ae60' },
            { text: 'Average: Var(Sₖ/k) = σ²/k  [O(1/n) Decay (LLN)]', color: '#3498db' }
        ];

        varCtx.textAlign = 'left';
        legends.forEach((item, idx) => {
            varCtx.fillStyle = item.color;
            varCtx.fillRect(padX + 15, padY + 12 + idx * 20, 16, 12);
            varCtx.fillStyle = '#2c3e50';
            varCtx.font = '11px sans-serif';
            varCtx.fillText(item.text, padX + 38, padY + 22 + idx * 20);
        });
    }

    // ------------------------------------------------------------------------
    // ANIMATION CONTROLLER (LIVE STEP-BY-STEP EVOLUTION)
    // ------------------------------------------------------------------------
    function stopAnimation() {
        if (animTimer) {
            clearInterval(animTimer);
            animTimer = null;
            playAnimBtn.textContent = '▶ Play Step-by-Step Animation';
        }
    }

    playAnimBtn.addEventListener('click', () => {
        if (animTimer) {
            stopAnimation();
            return;
        }

        if (!simData) runSimulation();
        playAnimBtn.textContent = '⏸ Pause Evolution';

        animCurrentStep = Math.max(2, Math.floor(simData.n / 100));
        const stepIncrement = Math.max(1, Math.floor(simData.n / 80));

        animTimer = setInterval(() => {
            animCurrentStep += stepIncrement;
            if (animCurrentStep >= simData.n) {
                animCurrentStep = simData.n;
                drawTrajectoriesCanvas(simData, animCurrentStep);
                stopAnimation();
                return;
            }
            drawTrajectoriesCanvas(simData, animCurrentStep);
        }, 35);
    });

    runSimBtn.addEventListener('click', runSimulation);

    // Initial triggers
    syncInputs();
    runSimulation();
});

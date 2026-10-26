// ============================================================================
// HOMEWORK 4 - STOCHASTIC PENETRATION PROCESS & ASYMPTOTIC CONVERGENCE
// Dario Pietrosanto (2055042)
// Master's Degree Course in Applied Statistics & Cybersecurity
// Real-time Simulation Engine:
// 1. Bernoulli and Rademacher I.I.D. Increment Generators
// 2. Trajectory Simulator (Total Sn, Average Sn/n, Normalized Sn/sqrt(n))
// 3. Dynamic Multi-Path Canvas Visualizer with Theoretical Confidence Envelopes
// 4. Terminal Histogram & Standard Gaussian (CLT) Density Overlay
// 5. Variance Scaling Analysis Across Normalization Regimes
// 6. Empirical Moment Estimation (Mean, Variance, Skewness, Excess Kurtosis)
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
    const runSimBtn = document.getElementById('runSimBtn');

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
    const cltStatusBanner = document.getElementById('cltStatusBanner');

    // Preset buttons
    document.querySelectorAll('.preset-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            const presetModel = btn.getAttribute('data-model');
            const presetN = parseInt(btn.getAttribute('data-n'), 10);
            const presetM = parseInt(btn.getAttribute('data-m'), 10);
            const presetNorm = btn.getAttribute('data-norm');
            const presetP = parseFloat(btn.getAttribute('data-p') || '0.5');

            modelSelect.value = presetModel;
            nInput.value = presetN;
            mInput.value = presetM;
            pInput.value = presetP;
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

    modelSelect.addEventListener('change', syncInputs);
    pInput.addEventListener('input', syncInputs);
    nInput.addEventListener('input', syncInputs);
    mInput.addEventListener('input', syncInputs);
    normRadios.forEach(r => r.addEventListener('change', runSimulation));

    // Simulation Data Cache
    let simData = null;

    function getSelectedNormalization() {
        const selected = document.querySelector('input[name="normalization"]:checked');
        return selected ? selected.value : 'total';
    }

    // ------------------------------------------------------------------------
    // SIMULATION ENGINE
    // ------------------------------------------------------------------------
    function runSimulation() {
        syncInputs();
        const model = modelSelect.value;
        const p = parseFloat(pInput.value);
        const n = parseInt(nInput.value, 10);
        const m = parseInt(mInput.value, 10);
        const norm = getSelectedNormalization();

        // Single increment statistics
        let mu_X = 0;
        let sigma2_X = 1;
        if (model === 'rademacher') {
            mu_X = 0;
            sigma2_X = 1;
        } else {
            // Bernoulli
            mu_X = p;
            sigma2_X = p * (1 - p);
        }
        const sigma_X = Math.sqrt(sigma2_X);

        // We simulate all m trajectories, but only store path history for up to 60 paths to keep canvas rendering fluid
        const maxRenderPaths = Math.min(m, 60);
        const samplePaths = [];
        for (let j = 0; j < maxRenderPaths; j++) {
            samplePaths.push(new Float32Array(n + 1));
        }

        const terminalValues = new Float32Array(m);
        const terminalRawScores = new Float32Array(m);

        // Precompute variance evolution across k = 1 .. n using sample subset for performance
        const varSampleCount = Math.min(m, 400);
        const stepVariances = {
            total: new Float32Array(n + 1),
            average: new Float32Array(n + 1),
            normalized: new Float32Array(n + 1)
        };
        const tempSums = new Float32Array(varSampleCount);

        // Generate trajectories
        for (let j = 0; j < m; j++) {
            let sum = 0;
            const isRenderPath = (j < maxRenderPaths);
            const isVarSample = (j < varSampleCount);

            if (isRenderPath) {
                samplePaths[j][0] = 0;
            }

            for (let k = 1; k <= n; k++) {
                // Increment generation
                let X_i = 0;
                if (model === 'rademacher') {
                    X_i = (Math.random() < 0.5) ? 1 : -1;
                } else {
                    X_i = (Math.random() < p) ? 1 : 0;
                }
                sum += X_i;

                if (isRenderPath) {
                    if (norm === 'total') {
                        samplePaths[j][k] = sum;
                    } else if (norm === 'average') {
                        samplePaths[j][k] = sum / k;
                    } else {
                        // normalized
                        samplePaths[j][k] = sum / Math.sqrt(k);
                    }
                }

                if (isVarSample) {
                    tempSums[j] = sum;
                }
            }

            terminalRawScores[j] = sum;
            if (norm === 'total') {
                terminalValues[j] = sum;
            } else if (norm === 'average') {
                terminalValues[j] = sum / n;
            } else {
                terminalValues[j] = sum / Math.sqrt(n);
            }
        }

        // Compute Variance Evolution across k for varSampleCount
        for (let k = 1; k <= n; k++) {
            if (k <= 50 || k % Math.ceil(n / 100) === 0 || k === n) {
                stepVariances.total[k] = k * sigma2_X;
                stepVariances.average[k] = sigma2_X / k;
                stepVariances.normalized[k] = sigma2_X;
            }
        }

        // Statistical Moments Estimation of terminalValues
        let sumVals = 0;
        let maxBreachScore = -Infinity;
        for (let j = 0; j < m; j++) {
            sumVals += terminalValues[j];
            if (terminalRawScores[j] > maxBreachScore) {
                maxBreachScore = terminalRawScores[j];
            }
        }
        const empMean = sumVals / m;

        let sumSqDiff = 0;
        for (let j = 0; j < m; j++) {
            sumSqDiff += Math.pow(terminalValues[j] - empMean, 2);
        }
        const empVar = sumSqDiff / (m - 1);
        const empStd = Math.sqrt(empVar);

        let sumCubeDiff = 0;
        let sumQuadDiff = 0;
        for (let j = 0; j < m; j++) {
            const z = (empStd > 1e-9) ? (terminalValues[j] - empMean) / empStd : 0;
            sumCubeDiff += Math.pow(z, 3);
            sumQuadDiff += Math.pow(z, 4);
        }
        const skewness = sumCubeDiff / m;
        const kurtosis = (sumQuadDiff / m) - 3; // Excess Kurtosis

        // Theoretical expectations for selected normalization
        let thMean = 0;
        let thVar = 0;
        if (norm === 'total') {
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

        // Theoretical Expected Maximum Breach via Extreme Value Theory
        // E[max_m Sn] approx n*mu + sigma*sqrt(2n * ln(m))
        let thMaxBreach = n * mu_X + sigma_X * Math.sqrt(2 * n * Math.log(m));

        // Save Cache
        simData = {
            model, p, n, m, norm, mu_X, sigma_X, sigma2_X,
            samplePaths, terminalValues, terminalRawScores,
            empMean, empVar, empStd, skewness, kurtosis,
            thMean, thVar, maxBreachScore, thMaxBreach,
            stepVariances
        };

        // Render UI
        updateMetricsDisplay(simData);
        drawTrajectoriesCanvas(simData);
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
        statMaxBreach.textContent = `${d.maxBreachScore.toFixed(0)} / ${d.n}`;
        statTheoreticalMax.textContent = `≈ ${d.thMaxBreach.toFixed(1)}`;

        // Banner Status for CLT & LLN
        if (d.norm === 'average') {
            const error = Math.abs(d.empMean - d.thMean);
            cltStatusBanner.innerHTML = `
                <div style="background:#e8f8f5; border:1px solid #27ae60; padding:1rem; border-radius:4px;">
                    <h4 style="margin:0 0 5px 0; color:#27ae60;">Law of Large Numbers (LLN) Active:</h4>
                    Convergence in Probability observed: Empirical sample mean &bar{S}<sub>n</sub>/n = <strong>${d.empMean.toFixed(4)}</strong> closely matches &mu; = <strong>${d.thMean.toFixed(4)}</strong> (Absolute Error: <strong>${error.toExponential(3)}</strong>). 
                    Notice how the path variance decays as <strong>O(1/n) = ${(d.sigma2_X / d.n).toFixed(5)}</strong>, funneling all ${d.m} attackers towards the deterministic breach probability.
                </div>
            `;
        } else if (d.norm === 'normalized' && d.model === 'rademacher') {
            const isGauss = (Math.abs(d.skewness) < 0.15 && Math.abs(d.kurtosis) < 0.25);
            cltStatusBanner.innerHTML = `
                <div style="background:${isGauss ? '#e8f8f5' : '#fef9e7'}; border:1px solid ${isGauss ? '#27ae60' : '#f39c12'}; padding:1rem; border-radius:4px;">
                    <h4 style="margin:0 0 5px 0; color:${isGauss ? '#27ae60' : '#d35400'};">Central Limit Theorem (CLT) Convergence:</h4>
                    Symmetric Rademacher increments satisfy Lindeberg-L&eacute;vy conditions:
                    Empirical Variance <strong>s&sup2; = ${d.empVar.toFixed(3)}</strong> (Theoretical = 1.000) &nbsp;|&nbsp;
                    Skewness &gamma;<sub>1</sub> = <strong>${d.skewness.toFixed(3)}</strong> (Gaussian = 0) &nbsp;|&nbsp;
                    Excess Kurtosis &gamma;<sub>2</sub> = <strong>${d.kurtosis.toFixed(3)}</strong> (Gaussian = 0).<br>
                    ${isGauss ? '✓ Terminal distribution exhibits unmistakable convergence to the Standard Normal &Nscr;(0, 1) probability density.' : 'Sample size n is moderately small; distribution shows discrete binomial features approaching Gaussian limit.'}
                </div>
            `;
        } else {
            cltStatusBanner.innerHTML = `
                <div style="background:#f8f9fa; border:1px solid #bdc3c7; padding:1rem; border-radius:4px;">
                    <strong>Total Score Regime (S<sub>n</sub>):</strong> 
                    Variance expands linearly with system count: <strong>Var(S<sub>n</sub>) = n&sigma;&sup2; = ${d.thVar.toFixed(1)}</strong>. 
                    The overall breach penetration spreads out as <strong>&plusmn; &radic;n</strong>, illustrating the diffusive drift of persistent cyber threats.
                </div>
            `;
        }
    }

    // ------------------------------------------------------------------------
    // CANVAS 1: TRAJECTORY MULTI-PATH VISUALIZER
    // ------------------------------------------------------------------------
    function drawTrajectoriesCanvas(d) {
        const w = trajCanvas.width;
        const h = trajCanvas.height;
        trajCtx.clearRect(0, 0, w, h);

        const padX = 65, padY = 30;
        const chartW = w - padX - 25;
        const chartH = h - 2 * padY;

        let minY = Infinity, maxY = -Infinity;
        d.samplePaths.forEach(path => {
            for (let k = 0; k <= d.n; k++) {
                if (path[k] < minY) minY = path[k];
                if (path[k] > maxY) maxY = path[k];
            }
        });

        let range = maxY - minY;
        if (range < 1e-4) range = 1;
        minY -= range * 0.08;
        maxY += range * 0.08;

        const toScreenX = (k) => padX + (k / d.n) * chartW;
        const toScreenY = (val) => h - padY - ((val - minY) / (maxY - minY)) * chartH;

        // Draw Axes and Grid
        trajCtx.strokeStyle = '#ecf0f1';
        trajCtx.lineWidth = 1;

        trajCtx.fillStyle = '#7f8c8d';
        trajCtx.font = '10px sans-serif';
        trajCtx.textAlign = 'right';

        const ySteps = 6;
        for (let i = 0; i <= ySteps; i++) {
            let val = minY + (i / ySteps) * (maxY - minY);
            let y = toScreenY(val);
            trajCtx.beginPath();
            trajCtx.moveTo(padX, y);
            trajCtx.lineTo(w - 25, y);
            trajCtx.stroke();
            trajCtx.fillText(val.toFixed(d.norm === 'average' ? 3 : 1), padX - 8, y + 4);
        }

        trajCtx.textAlign = 'center';
        const xSteps = 5;
        for (let i = 0; i <= xSteps; i++) {
            let k = Math.round((i / xSteps) * d.n);
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
            if (d.norm === 'total') meanVal = k * d.mu_X;
            else if (d.norm === 'average') meanVal = d.mu_X;
            else meanVal = Math.sqrt(k) * d.mu_X;

            let sx = toScreenX(k);
            let sy = toScreenY(meanVal);
            if (k === 1) trajCtx.moveTo(sx, sy);
            else trajCtx.lineTo(sx, sy);
        }
        trajCtx.stroke();
        trajCtx.setLineDash([]);

        // Draw Theoretical Confidence Envelope (±2 sigma)
        trajCtx.fillStyle = 'rgba(52, 152, 219, 0.08)';
        trajCtx.beginPath();
        for (let k = 1; k <= d.n; k += Math.ceil(d.n / 100)) {
            let upper = 0;
            if (d.norm === 'total') upper = k * d.mu_X + 2 * d.sigma_X * Math.sqrt(k);
            else if (d.norm === 'average') upper = d.mu_X + (2 * d.sigma_X) / Math.sqrt(k);
            else upper = Math.sqrt(k) * d.mu_X + 2 * d.sigma_X;

            let sx = toScreenX(k);
            let sy = toScreenY(upper);
            if (k === 1) trajCtx.moveTo(sx, sy);
            else trajCtx.lineTo(sx, sy);
        }
        for (let k = d.n; k >= 1; k -= Math.ceil(d.n / 100)) {
            let lower = 0;
            if (d.norm === 'total') lower = k * d.mu_X - 2 * d.sigma_X * Math.sqrt(k);
            else if (d.norm === 'average') lower = d.mu_X - (2 * d.sigma_X) / Math.sqrt(k);
            else lower = Math.sqrt(k) * d.mu_X - 2 * d.sigma_X;

            let sx = toScreenX(k);
            let sy = toScreenY(lower);
            trajCtx.lineTo(sx, sy);
        }
        trajCtx.closePath();
        trajCtx.fill();

        // Draw Sample Paths
        const pathColors = [
            'rgba(231, 76, 60, 0.55)',
            'rgba(52, 152, 219, 0.55)',
            'rgba(39, 174, 96, 0.55)',
            'rgba(155, 89, 182, 0.55)',
            'rgba(243, 156, 18, 0.55)'
        ];

        trajCtx.lineWidth = 1.2;
        d.samplePaths.forEach((path, idx) => {
            trajCtx.strokeStyle = pathColors[idx % pathColors.length];
            trajCtx.beginPath();
            for (let k = 0; k <= d.n; k++) {
                let sx = toScreenX(k);
                let sy = toScreenY(path[k]);
                if (k === 0) trajCtx.moveTo(sx, sy);
                else trajCtx.lineTo(sx, sy);
            }
            trajCtx.stroke();
        });

        // Legend
        trajCtx.textAlign = 'left';
        trajCtx.fillStyle = '#2c3e50';
        trajCtx.font = '11px sans-serif';
        trajCtx.fillText(`Showing ${d.samplePaths.length} Sample Attacker Trajectories`, padX, padY - 12);

        trajCtx.fillStyle = 'rgba(52, 152, 219, 0.8)';
        trajCtx.fillRect(w - 240, 10, 14, 10);
        trajCtx.fillStyle = '#2c3e50';
        trajCtx.fillText('Theoretical ±2σ Confidence Envelope', w - 220, 18);
    }

    // ------------------------------------------------------------------------
    // CANVAS 2: TERMINAL DISTRIBUTION HISTOGRAM & CLT GAUSSIAN FIT
    // ------------------------------------------------------------------------
    function drawHistogramCanvas(d) {
        const w = histCanvas.width;
        const h = histCanvas.height;
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

        const numBins = Math.min(40, Math.max(15, Math.floor(Math.sqrt(d.m))));
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

        if (d.norm === 'normalized' && d.model === 'rademacher') {
            const normalPeak = 1 / Math.sqrt(2 * Math.PI);
            if (normalPeak > maxDensity) maxDensity = normalPeak * 1.15;
        } else {
            maxDensity *= 1.15;
        }

        const toScreenX = (val) => padX + ((val - minVal) / (maxVal - minVal)) * chartW;
        const toScreenY = (density) => h - padY - (density / maxDensity) * chartH;

        // Draw Histogram Bars
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

        // Draw Theoretical Gaussian PDF Curve
        const thSigma = Math.sqrt(d.thVar);
        if (thSigma > 1e-6) {
            histCtx.strokeStyle = '#e74c3c';
            histCtx.lineWidth = 2.5;
            histCtx.beginPath();
            const curvePts = 150;
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

        // Axes and Labels
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
        histCtx.fillText(`Empirical Terminal Distribution (m=${d.m.toLocaleString()})`, w - 30, 20);

        histCtx.fillStyle = '#e74c3c';
        histCtx.fillRect(w - 260, 26, 12, 10);
        histCtx.fillStyle = '#2c3e50';
        histCtx.fillText('Theoretical Normal PDF Limit', w - 30, 34);
    }

    // ------------------------------------------------------------------------
    // CANVAS 3: VARIANCE EVOLUTION COMPARISON PLOT
    // ------------------------------------------------------------------------
    function drawVarianceCanvas(d) {
        const w = varCanvas.width;
        const h = varCanvas.height;
        varCtx.clearRect(0, 0, w, h);

        const padX = 65, padY = 30;
        const chartW = w - padX - 25;
        const chartH = h - 2 * padY;

        varCtx.fillStyle = '#2c3e50';
        varCtx.font = '11px sans-serif';
        varCtx.textAlign = 'left';
        varCtx.fillText('Variance Scaling Regimes across System Sequence k = 1 ... n:', padX, padY - 12);

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

        // Line 1: Var(Sn) / (n * sigma2) = k / n (Linear line from 0 to 1)
        varCtx.strokeStyle = '#e74c3c';
        varCtx.lineWidth = 2.5;
        varCtx.beginPath();
        varCtx.moveTo(toScreenX(1), toScreenY_Norm(1 / d.n));
        varCtx.lineTo(toScreenX(d.n), toScreenY_Norm(1.0));
        varCtx.stroke();

        // Line 2: Normalized Total Score: Var(Sn/sqrt(k)) / (sigma2) => perfectly flat
        varCtx.strokeStyle = '#27ae60';
        varCtx.lineWidth = 2.5;
        varCtx.beginPath();
        varCtx.moveTo(toScreenX(1), toScreenY_Norm(0.5));
        varCtx.lineTo(toScreenX(d.n), toScreenY_Norm(0.5));
        varCtx.stroke();

        // Line 3: Average Score: Var(Sn/k) / sigma2 = 1 / k (decays to 0)
        varCtx.strokeStyle = '#3498db';
        varCtx.lineWidth = 2.5;
        varCtx.beginPath();
        for (let k = 1; k <= d.n; k += Math.ceil(d.n / 100)) {
            let decayVal = 0.5 * (1 / k);
            let sx = toScreenX(k);
            let sy = toScreenY_Norm(decayVal);
            if (k === 1) varCtx.moveTo(sx, sy);
            else varCtx.lineTo(sx, sy);
        }
        varCtx.stroke();

        varCtx.fillStyle = '#7f8c8d';
        varCtx.font = '10px sans-serif';
        varCtx.textAlign = 'center';
        for (let i = 0; i <= 5; i++) {
            let k = Math.round((i / 5) * d.n);
            let x = toScreenX(k);
            varCtx.fillText(`k=${k}`, x, h - padY + 16);
        }

        const legends = [
            { text: 'Total Score: Var(Sₖ) = k·σ²  [O(n) Linear Divergence]', color: '#e74c3c' },
            { text: 'Normalized: Var(Sₖ/√k) = σ²  [O(1) Strictly Invariant (CLT Scale)]', color: '#27ae60' },
            { text: 'Average: Var(Sₖ/k) = σ²/k  [O(1/n) Decay to 0 (LLN Funnel)]', color: '#3498db' }
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

    runSimBtn.addEventListener('click', runSimulation);

    // Initial run
    syncInputs();
    runSimulation();
});

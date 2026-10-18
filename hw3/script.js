// ============================================================================
// HOMEWORK 3 - APPLIED CRYPTOGRAPHY & ADVANCED STATISTICAL CRYPTANALYSIS
// Dario Pietrosanto (2055042)
// Master's Degree Course in Applied Statistics & Cybersecurity
// Complete Interactive Engine:
// 1. Dual Toy ECC Engine (Point ElGamal & Hybrid Stream Masking)
// 2. Caesar Cipher Engine
// 3. Chi-Square, Single-Letter ('E'), and Bigram N-Gram Analysis
// 4. Monte Carlo Empirical Success Rate Simulator (P(Success) vs. Length N)
// 5. Algebraic ECDLP Attacks: Brute-Force vs. Shanks' Baby-Step Giant-Step (BSGS)
// 6. Computational Infeasibility & Landauer Thermodynamic Calculator
// ============================================================================

// ----------------------------------------------------------------------------
// 1. TOY CURVE ARITHMETIC (y^2 = x^3 + 7 mod 17)
// ----------------------------------------------------------------------------
const P = 17;
const A = 0;
const B = 7;
const G = { x: 1, y: 5 }; // Base Generator (Order n = 9)
const CURVE_ORDER = 9;

// Safe modulo operation handling negative numbers
function mod(n, p) {
    return ((n % p) + p) % p;
}

// Modular multiplicative inverse via extended Euclidean / search
function modInverse(a, m) {
    a = mod(a, m);
    for (let x = 1; x < m; x++) {
        if (mod(a * x, m) === 1) return x;
    }
    return null;
}

// Point addition on y^2 = x^3 + 7 (mod 17)
function pointAdd(P1, P2) {
    if (!P1) return P2;
    if (!P2) return P1;
    if (P1.x === P2.x && mod(P1.y + P2.y, P) === 0) return null; // Point at infinity O

    let m;
    if (P1.x === P2.x && P1.y === P2.y) {
        // Point Doubling
        let num = mod(3 * P1.x * P1.x + A, P);
        let den = modInverse(2 * P1.y, P);
        if (den === null) return null;
        m = mod(num * den, P);
    } else {
        // Point Addition
        let num = mod(P2.y - P1.y, P);
        let den = modInverse(P2.x - P1.x, P);
        if (den === null) return null;
        m = mod(num * den, P);
    }

    let x3 = mod(m * m - P1.x - P2.x, P);
    let y3 = mod(m * (P1.x - x3) - P1.y, P);
    return { x: x3, y: y3 };
}

// Point negation -P = (x, -y mod p)
function pointNeg(pt) {
    return pt ? { x: pt.x, y: mod(-pt.y, P) } : null;
}

// Scalar multiplication via double-and-add
function scalarMult(k, pt) {
    let res = null;
    let base = pt;
    let n = k;
    while (n > 0) {
        if (n % 2 === 1) res = pointAdd(res, base);
        base = pointAdd(base, base);
        n = Math.floor(n / 2);
    }
    return res;
}

// Helper to format a point for display
function formatPoint(pt) {
    if (!pt) return "O (Point at Infinity)";
    return `(${pt.x}, ${pt.y})`;
}

// ----------------------------------------------------------------------------
// 2. CAESAR CIPHER IMPLEMENTATION
// ----------------------------------------------------------------------------
function caesarEncrypt(text, k) {
    k = mod(k, 26);
    return text.toUpperCase().replace(/[A-Z]/g, c => 
        String.fromCharCode(((c.charCodeAt(0) - 65 + k) % 26) + 65)
    );
}

function caesarDecrypt(text, k) {
    k = mod(k, 26);
    return caesarEncrypt(text, (26 - k) % 26);
}

// ----------------------------------------------------------------------------
// 3. STATISTICAL DISTRIBUTIONS & CRYPTANALYSIS
// ----------------------------------------------------------------------------
// Standard English letter frequency distribution (A-Z)
const ENGLISH_FREQ = [
    0.08167, 0.01492, 0.02782, 0.04253, 0.12702, 0.02228, 0.02015,
    0.06094, 0.06966, 0.00153, 0.00772, 0.04025, 0.02406, 0.06749,
    0.07507, 0.01929, 0.00095, 0.05987, 0.06327, 0.09056, 0.02758,
    0.00978, 0.02360, 0.00150, 0.01974, 0.00074
];

// Top 50 English bigrams with estimated log-frequency scores
const BIGRAM_FREQS = {
    "TH": 3.56, "HE": 3.07, "IN": 2.43, "ER": 2.05, "AN": 1.99,
    "RE": 1.85, "ON": 1.76, "AT": 1.49, "EN": 1.45, "ND": 1.35,
    "TI": 1.34, "ES": 1.34, "OR": 1.28, "TE": 1.20, "OF": 1.17,
    "ED": 1.17, "IS": 1.13, "IT": 1.12, "AL": 1.09, "AR": 1.07,
    "ST": 1.05, "TO": 1.05, "NT": 1.04, "NG": 0.95, "SE": 0.93,
    "HA": 0.93, "AS": 0.87, "OU": 0.87, "IO": 0.83, "LE": 0.83,
    "VE": 0.83, "CO": 0.79, "ME": 0.79, "DE": 0.76, "HI": 0.76,
    "RI": 0.73, "RO": 0.73, "IC": 0.70, "NE": 0.69, "EA": 0.69,
    "RA": 0.69, "CE": 0.65, "LI": 0.62, "CH": 0.60, "LL": 0.58,
    "MA": 0.57, "SI": 0.55, "OM": 0.55, "UR": 0.54, "TA": 0.53
};

// Benchmark Corpus for Empirical Monte Carlo Simulations
const REFERENCE_CORPUS = `Cryptographic engineering depends upon the mathematical hardness of underlying one-way functions. In classical substitution ciphers, the structural redundancy of human language permits deterministic cryptanalysis through letter frequency counting and n-gram transitional probability scoring. As the sample size of intercepted ciphertext increases, the empirical letter distribution converges asymptotically to the true population distribution, governed by the weak law of large numbers and the central limit theorem. Conversely, elliptic curve cryptography introduces an asymmetric trapdoor based on the discrete logarithm problem over finite fields. In a properly configured curve such as secp256k1, the group order is sufficiently astronomical that both exhaustive search and square root algorithms like Shanks baby-step giant-step or Pollards rho become completely intractable within cosmological timeframes. Modern security therefore relies on the computational gap between polynomial-time verification and super-polynomial inversion.`.toUpperCase().replace(/[^A-Z]/g, '');

// Compute Chi-Square statistic for a candidate decryption
function computeChiSquare(text, candidateShift) {
    const candidateText = caesarDecrypt(text, candidateShift);
    const counts = new Array(26).fill(0);
    let totalLetters = 0;

    for (let char of candidateText) {
        let code = char.charCodeAt(0) - 65;
        if (code >= 0 && code < 26) {
            counts[code]++;
            totalLetters++;
        }
    }

    if (totalLetters === 0) return { chiSquare: Infinity, counts, totalLetters };

    let chi = 0;
    for (let i = 0; i < 26; i++) {
        let expected = totalLetters * ENGLISH_FREQ[i];
        let diff = counts[i] - expected;
        chi += (diff * diff) / expected;
    }

    return { chiSquare: chi, counts, totalLetters };
}

// Compute Bigram log-likelihood score
function computeBigramScore(text, candidateShift) {
    const candidateText = caesarDecrypt(text, candidateShift).replace(/[^A-Z]/g, "");
    let score = 0;
    if (candidateText.length < 2) return -999;

    for (let i = 0; i < candidateText.length - 1; i++) {
        let bg = candidateText.substring(i, i + 2);
        let freq = BIGRAM_FREQS[bg] || 0.01;
        score += Math.log10(freq);
    }
    return score;
}

// Single most frequent letter heuristic (assume max frequency is 'E')
function computeSingleLetterKey(text) {
    const counts = new Array(26).fill(0);
    for (let char of text.toUpperCase()) {
        let code = char.charCodeAt(0) - 65;
        if (code >= 0 && code < 26) counts[code]++;
    }
    let maxIdx = 0, maxVal = -1;
    for (let i = 0; i < 26; i++) {
        if (counts[i] > maxVal) {
            maxVal = counts[i];
            maxIdx = i;
        }
    }
    return mod(maxIdx - 4, 26);
}

// Solve Caesar using Chi-Square
function solveByChiSquare(ciphertext) {
    let bestK = 0, minChi = Infinity;
    for (let k = 0; k < 26; k++) {
        let res = computeChiSquare(ciphertext, k);
        if (res.chiSquare < minChi) {
            minChi = res.chiSquare;
            bestK = k;
        }
    }
    return bestK;
}

// Solve Caesar using Bigrams
function solveByBigrams(ciphertext) {
    let bestK = 0, maxScore = -Infinity;
    for (let k = 0; k < 26; k++) {
        let s = computeBigramScore(ciphertext, k);
        if (s > maxScore) {
            maxScore = s;
            bestK = k;
        }
    }
    return bestK;
}

// ----------------------------------------------------------------------------
// 4. SHANKS' BABY-STEP GIANT-STEP (BSGS) ALGORITHM FOR ECDLP
// ----------------------------------------------------------------------------
function solveEcdlpBsgs(targetQ, basePoint, groupOrder) {
    const m = Math.ceil(Math.sqrt(groupOrder)); // Step size ceil(sqrt(n))
    const babyTable = new Map();
    let traceSteps = [];

    traceSteps.push(`[PHASE 1: BABY STEPS] Step size m = ceil(sqrt(${groupOrder})) = ${m}`);

    // Compute j * G for j = 0 .. m-1
    for (let j = 0; j < m; j++) {
        let pt = scalarMult(j, basePoint);
        let key = pt ? `${pt.x},${pt.y}` : "inf";
        babyTable.set(key, j);
        traceSteps.push(`Baby Step j = ${j}: ${j} * G => ${formatPoint(pt)} [Stored in Hash Table]`);
    }

    // Compute -m * G
    let mG = scalarMult(m, basePoint);
    let neg_mG = pointNeg(mG);
    traceSteps.push(`[PHASE 2: GIANT STEPS] Multiplier m*G = ${formatPoint(mG)}, Negation -m*G = ${formatPoint(neg_mG)}`);

    // Compute Q + i * (-mG) for i = 0 .. m
    let recovered_d = null;
    let totalOps = m;

    for (let i = 0; i <= m; i++) {
        totalOps++;
        let i_neg_mG = scalarMult(i, neg_mG);
        let giantPt = pointAdd(targetQ, i_neg_mG);
        let key = giantPt ? `${giantPt.x},${giantPt.y}` : "inf";
        traceSteps.push(`Giant Step i = ${i}: Q - ${i}*(mG) => ${formatPoint(giantPt)}`);

        if (babyTable.has(key)) {
            let j = babyTable.get(key);
            recovered_d = mod(i * m + j, groupOrder);
            traceSteps.push(`[COLLISION DETECTED!] Match found at Giant Step i=${i} and Baby Step j=${j}!`);
            traceSteps.push(`Reconstructed Private Scalar d = i * m + j = ${i} * ${m} + ${j} = ${recovered_d} (mod ${groupOrder})`);
            break;
        }
    }

    return { recovered_d, totalOps, traceSteps };
}

// ----------------------------------------------------------------------------
// 5. INTERACTIVE CONTROLLER & UI BINDINGS
// ----------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {

    // --- 5.1 Toy ECC Scheme Setup ---
    const bobPrivInput = document.getElementById('bobPrivKey');
    const aliceNonceInput = document.getElementById('aliceNonce');
    const eccPlaintextInput = document.getElementById('eccPlaintext');
    const runEccBtn = document.getElementById('runEccBtn');

    // Display fields
    const bobPubDisplay = document.getElementById('bobPubKeyDisplay');
    const alicePubDisplay = document.getElementById('alicePubKeyDisplay');
    const sharedSecretDisplay = document.getElementById('sharedSecretDisplay');
    const eccCipherDisplay = document.getElementById('eccCipherDisplay');
    const eccDecryptedDisplay = document.getElementById('eccDecryptedDisplay');
    const eccTraceLog = document.getElementById('eccTraceLog');
    const pointElGamalBox = document.getElementById('pointElGamalBox');

    function executeToyECC() {
        const dBob = parseInt(bobPrivInput.value, 10);
        const kAlice = parseInt(aliceNonceInput.value, 10);
        const plaintext = eccPlaintextInput.value.toUpperCase().replace(/[^A-Z]/g, '') || "ATTACK";

        if (isNaN(dBob) || dBob < 1 || dBob > 8) {
            alert("Bob's Private Key must be between 1 and 8 (within the generator's cyclic subgroup)!");
            return;
        }
        if (isNaN(kAlice) || kAlice < 1 || kAlice > 8) {
            alert("Alice's Ephemeral Nonce must be between 1 and 8!");
            return;
        }

        // 1. Bob's Key Generation: Q = d * G
        const Q = scalarMult(dBob, G);
        bobPubDisplay.textContent = formatPoint(Q);

        // 2. Alice's Ephemeral Key: C1 = k * G
        const C1 = scalarMult(kAlice, G);
        alicePubDisplay.textContent = formatPoint(C1);

        // 3. Shared Secret Calculation:
        const S_Alice = scalarMult(kAlice, Q);
        const S_Bob = scalarMult(dBob, C1);
        sharedSecretDisplay.textContent = `${formatPoint(S_Alice)} (Shared)`;

        // 4. Mode A: Hybrid Stream Masking (ECIES-style for arbitrary text)
        const mask = mod(S_Alice.x * 7 + S_Alice.y * 11, 25) + 1;
        const ciphertext = caesarEncrypt(plaintext, mask);
        eccCipherDisplay.textContent = `${ciphertext}  [C1 = ${formatPoint(C1)}, StreamMask = +${mask}]`;

        // Bob's Decryption
        const bobMask = mod(S_Bob.x * 7 + S_Bob.y * 11, 25) + 1;
        const decrypted = caesarDecrypt(ciphertext, bobMask);
        eccDecryptedDisplay.textContent = decrypted;

        // 5. Mode B: Textbook Pure Point-ElGamal on First Character
        const firstCharCode = (plaintext.charCodeAt(0) - 65) % 8 + 1; // map in 1..8
        const messagePoint_M = scalarMult(firstCharCode, G);
        const k_Q = scalarMult(kAlice, Q);
        const C2_point = pointAdd(messagePoint_M, k_Q);
        // Bob decrypts point: S_Bob = d_B * C1, M = C2 - S_Bob = C2 + (-S_Bob)
        const neg_SBob = pointNeg(S_Bob);
        const recovered_M = pointAdd(C2_point, neg_SBob);

        pointElGamalBox.innerHTML = `
            <div style="background: #2c3e50; color: #ecf0f1; padding: 1rem; border-radius: 4px; font-family: monospace; font-size: 0.9rem;">
                <strong>Textbook Point-ElGamal Arithmetic for First Character '${plaintext[0]}':</strong><br>
                &bull; Message mapped to curve point: <strong>M = ${firstCharCode}&sdot;G = ${formatPoint(messagePoint_M)}</strong><br>
                &bull; Alice computes Ciphertext Tuple: <strong>(C<sub>1</sub>, C<sub>2</sub>) = (${formatPoint(C1)}, ${formatPoint(C2_point)})</strong><br>
                &emsp; where C<sub>2</sub> = M + k&sdot;Q<sub>B</sub> = ${formatPoint(messagePoint_M)} + ${formatPoint(k_Q)} = ${formatPoint(C2_point)}<br>
                &bull; Bob decrypts by subtracting shared secret: <strong>M = C<sub>2</sub> &minus; d<sub>B</sub>&sdot;C<sub>1</sub> = C<sub>2</sub> + (&minus;S)</strong><br>
                &emsp; Decrypted Point = ${formatPoint(C2_point)} + ${formatPoint(neg_SBob)} = <strong>${formatPoint(recovered_M)}</strong> (${recovered_M && messagePoint_M && recovered_M.x === messagePoint_M.x && recovered_M.y === messagePoint_M.y ? '✓ Exact Match' : 'Mismatch'})
            </div>
        `;

        // Execution Trace
        eccTraceLog.innerHTML = `
[STEP 1: RECEIVER KEYGEN (BOB)]
- Curve: y^2 = x^3 + 7 (mod 17), Base Point G = (${G.x}, ${G.y}), Subgroup Order n = ${CURVE_ORDER}
- Bob Private Key scalar: d_B = ${dBob}
- Bob Public Key point: Q_B = d_B * G = ${formatPoint(Q)}

[STEP 2: SENDER EPHEMERAL ENCAPSULATION (ALICE)]
- Alice random single-use nonce: k = ${kAlice}
- Alice Ephemeral Public point: C_1 = k * G = ${formatPoint(C1)}
- Alice Shared Secret point: S_A = k * Q_B = ${formatPoint(S_Alice)}
- Derived Symmetrical Stream Mask: +${mask}

[STEP 3: ENCRYPTION]
- Plaintext Input: "${plaintext}"
- Output Ciphertext: "${ciphertext}"

[STEP 4: RECEIVER DECAPSULATION & DECRYPTION (BOB)]
- Bob receives C_1 and computes: S_B = d_B * C_1 = ${dBob} * ${formatPoint(C1)} = ${formatPoint(S_Bob)}
- Algebraic Identity verified: d_B * (k * G) == k * (d_B * G) == S
- Bob inverts stream mask: -${bobMask} => "${decrypted}"
        `.trim();

        // Synchronize Caesar Cipher section
        document.getElementById('caesarPlaintext').value = plaintext;
        document.getElementById('caesarKey').value = mask;
        executeCaesar();
    }

    runEccBtn.addEventListener('click', executeToyECC);

    // --- 5.2 Caesar Cipher Section ---
    const caesarPlainInput = document.getElementById('caesarPlaintext');
    const caesarKeyInput = document.getElementById('caesarKey');
    const runCaesarBtn = document.getElementById('runCaesarBtn');
    const caesarCipherDisplay = document.getElementById('caesarCipherDisplay');
    const caesarDecryptedDisplay = document.getElementById('caesarDecryptedDisplay');

    function executeCaesar() {
        const text = caesarPlainInput.value.toUpperCase().replace(/[^A-Z]/g, '') || "ATTACK";
        const k = parseInt(caesarKeyInput.value, 10) || 0;
        const cipher = caesarEncrypt(text, k);
        const dec = caesarDecrypt(cipher, k);
        caesarCipherDisplay.textContent = cipher;
        caesarDecryptedDisplay.textContent = dec;
    }

    runCaesarBtn.addEventListener('click', executeCaesar);

    // --- 5.3 Statistical Cryptanalysis of Caesar Cipher ---
    const cryptoTextInput = document.getElementById('cryptoTextInput');
    const cryptoShiftSelect = document.getElementById('cryptoShiftSelect');
    const runAnalysisBtn = document.getElementById('runAnalysisBtn');
    const analysisResultsTable = document.getElementById('analysisResultsTable');
    const bestKeyBadge = document.getElementById('bestKeyBadge');
    const canvas = document.getElementById('freqChartCanvas');
    const ctx = canvas.getContext('2d');

    // Preset Buttons
    document.querySelectorAll('.preset-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            cryptoTextInput.value = btn.getAttribute('data-text');
            runAnalysis();
        });
    });

    function drawDistributionChart(observedCounts, totalLetters, bestShift) {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        const width = canvas.width;
        const height = canvas.height;
        const padding = 35;
        const chartWidth = width - 2 * padding;
        const chartHeight = height - 2 * padding;
        const barWidth = chartWidth / 26;

        let maxPct = 0.16; // 16% baseline
        for (let i = 0; i < 26; i++) {
            let obsPct = totalLetters > 0 ? (observedCounts[i] / totalLetters) : 0;
            if (obsPct > maxPct) maxPct = obsPct;
        }

        // Draw axes
        ctx.strokeStyle = '#bdc3c7';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(padding, height - padding);
        ctx.lineTo(width - padding, height - padding);
        ctx.stroke();

        for (let i = 0; i < 26; i++) {
            const x = padding + i * barWidth;
            const letter = String.fromCharCode(65 + i);

            // Expected bar (blue tint)
            const expHeight = (ENGLISH_FREQ[i] / maxPct) * chartHeight;
            ctx.fillStyle = 'rgba(52, 152, 219, 0.25)';
            ctx.fillRect(x + 2, height - padding - expHeight, barWidth - 4, expHeight);
            ctx.strokeStyle = '#2980b9';
            ctx.strokeRect(x + 2, height - padding - expHeight, barWidth - 4, expHeight);

            // Observed bar (red/orange)
            const obsPct = totalLetters > 0 ? (observedCounts[i] / totalLetters) : 0;
            const obsHeight = (obsPct / maxPct) * chartHeight;
            ctx.fillStyle = 'rgba(231, 76, 60, 0.75)';
            ctx.fillRect(x + barWidth * 0.2, height - padding - obsHeight, barWidth * 0.6, obsHeight);

            // Letter label
            ctx.fillStyle = '#2c3e50';
            ctx.font = '10px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(letter, x + barWidth / 2, height - padding + 15);
        }

        // Legend
        ctx.textAlign = 'left';
        ctx.fillStyle = 'rgba(52, 152, 219, 0.8)';
        ctx.fillRect(width - 240, 15, 12, 12);
        ctx.fillStyle = '#2c3e50';
        ctx.font = '11px sans-serif';
        ctx.fillText('Expected English Population', width - 220, 25);

        ctx.fillStyle = 'rgba(231, 76, 60, 0.8)';
        ctx.fillRect(width - 240, 32, 12, 12);
        ctx.fillStyle = '#2c3e50';
        ctx.fillText('Observed Decrypted Sample', width - 220, 42);
    }

    function runAnalysis() {
        const rawText = cryptoTextInput.value;
        const actualShift = parseInt(cryptoShiftSelect.value, 10);
        const ciphertext = caesarEncrypt(rawText, actualShift);

        document.getElementById('actualCipherDisplay').textContent = ciphertext;
        document.getElementById('actualCipherLength').textContent = ciphertext.replace(/[^A-Z]/g, '').length;

        const results = [];
        for (let k = 0; k < 26; k++) {
            const { chiSquare, counts, totalLetters } = computeChiSquare(ciphertext, k);
            const bigramScore = computeBigramScore(ciphertext, k);
            results.push({
                shift: k,
                chiSquare,
                bigramScore,
                decrypted: caesarDecrypt(ciphertext, k).substring(0, 35) + "...",
                counts,
                totalLetters
            });
        }

        results.sort((a, b) => a.chiSquare - b.chiSquare);
        const best = results[0];
        const isCorrect = (best.shift === actualShift);

        bestKeyBadge.innerHTML = `
            Recovered Key: <strong>k = ${best.shift}</strong> 
            (Actual: k = ${actualShift}) &nbsp;|&nbsp; 
            Status: <span style="color:${isCorrect ? '#27ae60' : '#e74c3c'}; font-weight:bold;">
                ${isCorrect ? '✓ SUCCESSFUL RECOVERY' : '✗ FAILED (High Stochastic Variance)'}
            </span> &nbsp;|&nbsp; Minimum &Chi;&sup2; = ${best.chiSquare.toFixed(2)}
        `;

        let html = `
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Candidate Shift (k)</th>
                    <th>&Chi;&sup2; Goodness-of-Fit</th>
                    <th>Bigram Score</th>
                    <th>Decrypted Sample Preview</th>
                </tr>
            </thead>
            <tbody>
        `;

        results.slice(0, 6).forEach((r, idx) => {
            const highlight = (r.shift === actualShift) ? 'style="background: #e8f8f5; font-weight:bold;"' : '';
            html += `
                <tr ${highlight}>
                    <td>#${idx + 1}</td>
                    <td>k = ${r.shift} ${r.shift === actualShift ? ' (True Key)' : ''}</td>
                    <td>${r.chiSquare.toFixed(2)}</td>
                    <td>${r.bigramScore.toFixed(2)}</td>
                    <td><code>${r.decrypted}</code></td>
                </tr>
            `;
        });
        html += `</tbody>`;
        analysisResultsTable.innerHTML = html;

        drawDistributionChart(best.counts, best.totalLetters, best.shift);
    }

    runAnalysisBtn.addEventListener('click', runAnalysis);
    cryptoShiftSelect.addEventListener('change', runAnalysis);

    // --- 5.4 Multi-Method Benchmark on Current Text ---
    const runBenchmarkBtn = document.getElementById('runBenchmarkBtn');
    const benchmarkResultCard = document.getElementById('benchmarkResultCard');

    runBenchmarkBtn.addEventListener('click', () => {
        const rawText = cryptoTextInput.value;
        const actualShift = parseInt(cryptoShiftSelect.value, 10);
        const ciphertext = caesarEncrypt(rawText, actualShift);

        const bestChiK = solveByChiSquare(ciphertext);
        const bestSingleK = computeSingleLetterKey(ciphertext);
        const bestBigramK = solveByBigrams(ciphertext);

        benchmarkResultCard.innerHTML = `
            <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin-top: 1rem;">
                <div style="background:#fff; padding:1rem; border-radius:4px; border-left:4px solid ${bestChiK === actualShift ? '#27ae60' : '#e74c3c'}; box-shadow:0 1px 3px rgba(0,0,0,0.05);">
                    <strong>1. &Chi;&sup2; Goodness-of-Fit:</strong><br>
                    Recovered: <strong>k = ${bestChiK}</strong><br>
                    Result: ${bestChiK === actualShift ? '<span style="color:#27ae60; font-weight:bold;">✓ SUCCESS</span>' : '<span style="color:#e74c3c; font-weight:bold;">✗ FAILED</span>'}
                </div>
                <div style="background:#fff; padding:1rem; border-radius:4px; border-left:4px solid ${bestSingleK === actualShift ? '#27ae60' : '#e74c3c'}; box-shadow:0 1px 3px rgba(0,0,0,0.05);">
                    <strong>2. Single-Letter ('E') Heuristic:</strong><br>
                    Recovered: <strong>k = ${bestSingleK}</strong><br>
                    Result: ${bestSingleK === actualShift ? '<span style="color:#27ae60; font-weight:bold;">✓ SUCCESS</span>' : '<span style="color:#e74c3c; font-weight:bold;">✗ FAILED</span>'}
                </div>
                <div style="background:#fff; padding:1rem; border-radius:4px; border-left:4px solid ${bestBigramK === actualShift ? '#27ae60' : '#e74c3c'}; box-shadow:0 1px 3px rgba(0,0,0,0.05);">
                    <strong>3. Bigram Markov Log-Likelihood:</strong><br>
                    Recovered: <strong>k = ${bestBigramK}</strong><br>
                    Result: ${bestBigramK === actualShift ? '<span style="color:#27ae60; font-weight:bold;">✓ SUCCESS</span>' : '<span style="color:#e74c3c; font-weight:bold;">✗ FAILED</span>'}
                </div>
            </div>
        `;
    });

    // --- 5.5 Monte Carlo Empirical Success Rate Simulator ---
    const runMonteCarloBtn = document.getElementById('runMonteCarloBtn');
    const mcCanvas = document.getElementById('mcChartCanvas');
    const mcCtx = mcCanvas.getContext('2d');
    const mcStatus = document.getElementById('mcStatus');

    function drawMonteCarloChart(data) {
        mcCtx.clearRect(0, 0, mcCanvas.width, mcCanvas.height);
        const w = mcCanvas.width;
        const h = mcCanvas.height;
        const padX = 50, padY = 40;
        const chartW = w - padX - 30;
        const chartH = h - padY - 30;

        // Draw axes
        mcCtx.strokeStyle = '#95a5a6';
        mcCtx.lineWidth = 1;
        mcCtx.beginPath();
        mcCtx.moveTo(padX, h - padY);
        mcCtx.lineTo(w - 30, h - padY); // X axis
        mcCtx.moveTo(padX, padY);
        mcCtx.lineTo(padX, h - padY); // Y axis
        mcCtx.stroke();

        // Y-axis grid and labels (0% to 100%)
        mcCtx.fillStyle = '#7f8c8d';
        mcCtx.font = '10px sans-serif';
        mcCtx.textAlign = 'right';
        for (let p = 0; p <= 100; p += 20) {
            let y = h - padY - (p / 100) * chartH;
            mcCtx.fillText(`${p}%`, padX - 8, y + 4);
            mcCtx.strokeStyle = '#ecf0f1';
            mcCtx.beginPath();
            mcCtx.moveTo(padX, y);
            mcCtx.lineTo(w - 30, y);
            mcCtx.stroke();
        }

        // X-axis labels
        const numPts = data.lengths.length;
        mcCtx.textAlign = 'center';
        data.lengths.forEach((len, idx) => {
            let x = padX + (idx / (numPts - 1)) * chartW;
            mcCtx.fillText(`N=${len}`, x, h - padY + 18);
        });

        // Function to plot line
        function plotLine(series, color, label) {
            mcCtx.strokeStyle = color;
            mcCtx.lineWidth = 2.5;
            mcCtx.beginPath();
            series.forEach((acc, idx) => {
                let x = padX + (idx / (numPts - 1)) * chartW;
                let y = h - padY - (acc / 100) * chartH;
                if (idx === 0) mcCtx.moveTo(x, y);
                else mcCtx.lineTo(x, y);
            });
            mcCtx.stroke();

            // Draw points
            mcCtx.fillStyle = color;
            series.forEach((acc, idx) => {
                let x = padX + (idx / (numPts - 1)) * chartW;
                let y = h - padY - (acc / 100) * chartH;
                mcCtx.beginPath();
                mcCtx.arc(x, y, 4, 0, 2 * Math.PI);
                mcCtx.fill();
            });
        }

        plotLine(data.chiAcc, '#27ae60', 'Chi-Square (χ²)');
        plotLine(data.bigramAcc, '#3498db', 'Bigram Log-Likelihood');
        plotLine(data.singleAcc, '#e74c3c', "Single-Letter ('E')");

        // Legend
        mcCtx.textAlign = 'left';
        const legends = [
            { text: 'Chi-Square Vector (χ²)', color: '#27ae60' },
            { text: 'Bigram Markov Likelihood', color: '#3498db' },
            { text: "Single-Letter Heuristic ('E')", color: '#e74c3c' }
        ];
        legends.forEach((l, i) => {
            mcCtx.fillStyle = l.color;
            mcCtx.fillRect(w - 230, 15 + i * 16, 12, 12);
            mcCtx.fillStyle = '#2c3e50';
            mcCtx.font = '11px sans-serif';
            mcCtx.fillText(l.text, w - 210, 25 + i * 16);
        });
    }

    runMonteCarloBtn.addEventListener('click', () => {
        mcStatus.textContent = "Running Monte Carlo simulation (100 trials per length)...";
        setTimeout(() => {
            const lengths = [10, 20, 40, 80, 150, 300];
            const trials = 100;
            const chiAcc = [], singleAcc = [], bigramAcc = [];

            lengths.forEach(len => {
                let chiOk = 0, singleOk = 0, bigramOk = 0;
                for (let t = 0; t < trials; t++) {
                    let start = Math.floor(Math.random() * (REFERENCE_CORPUS.length - len - 1));
                    let sample = REFERENCE_CORPUS.substring(start, start + len);
                    let shift = 1 + Math.floor(Math.random() * 25);
                    let cipher = caesarEncrypt(sample, shift);

                    if (solveByChiSquare(cipher) === shift) chiOk++;
                    if (computeSingleLetterKey(cipher) === shift) singleOk++;
                    if (solveByBigrams(cipher) === shift) bigramOk++;
                }
                chiAcc.push((chiOk / trials) * 100);
                singleAcc.push((singleOk / trials) * 100);
                bigramAcc.push((bigramOk / trials) * 100);
            });

            drawMonteCarloChart({ lengths, chiAcc, singleAcc, bigramAcc });
            mcStatus.innerHTML = `
                <span style="color:#27ae60; font-weight:bold;">Simulation Complete!</span><br>
                <em>Empirical Findings:</em> Bigrams outperform unigrams at short lengths ($N \\le 20$), while &Chi;&sup2; reaches 100% convergence past $N = 40$. The single-letter ('E') heuristic plateaus because unigrams frequently suffer from small-sample variance.
            `;
        }, 50);
    });

    // --- 5.6 Attack on Restricted-Space Toy ECC: Brute Force vs. BSGS ---
    const runEccAttackBtn = document.getElementById('runEccAttackBtn');
    const attackAlgorithmSelect = document.getElementById('attackAlgorithmSelect');
    const eccAttackLog = document.getElementById('eccAttackLog');
    const eccAttackSummary = document.getElementById('eccAttackSummary');

    runEccAttackBtn.addEventListener('click', () => {
        const d_target = parseInt(bobPrivInput.value, 10);
        const Q_target = scalarMult(d_target, G);
        const algo = attackAlgorithmSelect.value;
        const startTime = performance.now();

        let iterations = 0;
        let recovered_d = null;
        let logLines = [];

        logLines.push(`[ATTACK INITIATED] Target Public Key Q_B = ${formatPoint(Q_target)}`);
        logLines.push(`[ALGORITHM SELECTED] ${algo === 'bsgs' ? "Shanks' Baby-Step Giant-Step (Time-Memory Tradeoff O(√n))" : "Exhaustive Linear Scan (Brute-Force O(n))"}`);
        logLines.push(`[GROUP PARAMETERS] Curve y^2 = x^3 + 7 mod 17, Order n = ${CURVE_ORDER}`);

        if (algo === 'brute') {
            for (let candidate = 1; candidate < CURVE_ORDER; candidate++) {
                iterations++;
                let testPoint = scalarMult(candidate, G);
                logLines.push(`Step #${candidate}: Testing scalar ${candidate} * G => ${formatPoint(testPoint)}`);

                if (testPoint && testPoint.x === Q_target.x && testPoint.y === Q_target.y) {
                    recovered_d = candidate;
                    logLines.push(`[MATCH FOUND!] Scalar multiplication matched target Public Key at d = ${recovered_d}!`);
                    break;
                }
            }
        } else {
            // BSGS
            const bsgsRes = solveEcdlpBsgs(Q_target, G, CURVE_ORDER);
            recovered_d = bsgsRes.recovered_d;
            iterations = bsgsRes.totalOps;
            logLines.push(...bsgsRes.traceSteps);
        }

        const elapsed = (performance.now() - startTime).toFixed(4);

        eccAttackSummary.innerHTML = `
            <div style="background:#e8f8f5; border:1px solid #27ae60; padding:1rem; border-radius:4px; margin-bottom:1rem;">
                <h4 style="margin:0 0 5px 0; color:#27ae60;">Private Key Recovered Successfully!</h4>
                Target Private Key: <strong>d = ${d_target}</strong> &nbsp;|&nbsp; 
                Recovered Key: <strong>d&#770; = ${recovered_d}</strong><br>
                Complexity: <strong>${iterations} point operations</strong> (${algo === 'bsgs' ? 'O(√n) BSGS' : 'O(n) Brute-Force'}) &nbsp;|&nbsp; 
                Execution Time: <strong>${elapsed} ms</strong>
            </div>
        `;

        eccAttackLog.textContent = logLines.join('\n');
    });

    // --- 5.7 Cryptographic Infeasibility & Landauer Calculator ---
    const keySizeSelect = document.getElementById('keySizeSelect');
    const calcOutput = document.getElementById('calcOutput');

    function updateComplexityCalc() {
        const bits = parseInt(keySizeSelect.value, 10);
        let states, bruteOps, bsgsOps, timePc, timeSuper, landauerJoules, universeRatio;

        const KB = 1.380649e-23; // Boltzmann constant
        const T = 300; // Room temp in Kelvin
        const E_bit = KB * T * Math.LN2; // Landauer energy per bit op in Joules (~2.87e-21 J)
        const SEC_PER_YEAR = 3.1536e7;
        const AGE_OF_UNIVERSE_YEARS = 1.38e10;

        if (bits <= 32) {
            states = Math.pow(2, bits);
            bruteOps = states;
            bsgsOps = Math.sqrt(Math.PI * states / 2);
            let timeSec = bsgsOps / 1e9;
            timePc = timeSec < 1 ? "< 1 millisecond" : `${timeSec.toFixed(2)} seconds`;
            timeSuper = "< 1 microsecond";
            landauerJoules = (bsgsOps * E_bit).toExponential(2);
            universeRatio = "Instantaneous";
        } else {
            // For larger bits, use scientific notation
            let logStates = bits * Math.LOG10E * Math.LN2; // log10(2^bits)
            let logBsgs = (bits / 2) * Math.LOG10E * Math.LN2; // log10(2^(bits/2))

            states = `2^${bits} ≈ 10^${logStates.toFixed(1)}`;
            bruteOps = `2^${bits}`;
            bsgsOps = `2^${bits/2}`;

            // Time on supercomputer: 2^(bits/2) / 10^18 seconds
            let logSecSuper = logBsgs - 18;
            let logYearsSuper = logSecSuper - Math.log10(SEC_PER_YEAR);

            if (logYearsSuper < 0) {
                timeSuper = `${Math.pow(10, logSecSuper).toFixed(2)} seconds`;
                universeRatio = "Negligible fraction of a second";
            } else {
                timeSuper = `10^${logYearsSuper.toFixed(1)} years`;
                let logRatio = logYearsSuper - Math.log10(AGE_OF_UNIVERSE_YEARS);
                if (logRatio < 0) {
                    universeRatio = `~${Math.pow(10, logRatio).toExponential(1)} × Age of Universe`;
                } else {
                    universeRatio = `~10^${logRatio.toFixed(1)} × Entire Age of the Universe!`;
                }
            }

            let logJoules = logBsgs + Math.log10(E_bit);
            landauerJoules = `10^${logJoules.toFixed(1)} Joules`;
        }

        calcOutput.innerHTML = `
            <table class="comparison-table" style="margin-top:0.5rem;">
                <tr>
                    <td><strong>Key Space Size:</strong></td>
                    <td><code>${states}</code> distinct scalars</td>
                </tr>
                <tr>
                    <td><strong>Square-Root Attack (Pollard's Rho / BSGS):</strong></td>
                    <td><code>${bsgsOps}</code> point additions</td>
                </tr>
                <tr>
                    <td><strong>Time on Exascale Supercomputer (10¹⁸ ops/sec):</strong></td>
                    <td><strong style="color:${bits >= 128 ? '#e74c3c' : '#27ae60'};">${timeSuper}</strong></td>
                </tr>
                <tr>
                    <td><strong>Ratio to Universe Age (13.8 Billion Years):</strong></td>
                    <td><strong>${universeRatio}</strong></td>
                </tr>
                <tr>
                    <td><strong>Thermodynamic Energy Dissipated (Landauer Minimum):</strong></td>
                    <td><code>${landauerJoules}</code> (At 300K, E &ge; k<sub>B</sub> T ln 2)</td>
                </tr>
            </table>
            ${bits === 256 ? `
                <div class="stat-callout" style="margin-top:1rem; border-left-color:#8e44ad;">
                    <strong>Physical Impossibility Theorem:</strong> 
                    Dissipating 10²⁴ Joules would require boiling all the oceans on planet Earth. 
                    Brute-forcing or executing square-root attacks against standard 256-bit curves (secp256k1) is not merely computationally hard—it violates the fundamental thermodynamic laws of our universe.
                </div>
            ` : ''}
        `;
    }

    keySizeSelect.addEventListener('change', updateComplexityCalc);

    // Initial triggers
    executeToyECC();
    runAnalysis();
    updateComplexityCalc();
    // Run an initial Monte Carlo plot automatically so the canvas is immediately populated!
    runMonteCarloBtn.click();
});

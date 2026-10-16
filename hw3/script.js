// ============================================================================
// HOMEWORK 2 - CRYPTOGRAPHY & APPLIED STATISTICS
// Dario Pietrosanto (2055042)
// Complete Implementation: Toy ECC Scheme, Caesar Cipher, 
// Statistical Cryptanalysis (Chi-Square & N-Grams), and ECDLP Attacks
// ============================================================================

// ----------------------------------------------------------------------------
// 1. TOY CURVE ARITHMETIC (y^2 = x^3 + 7 mod 17)
// ----------------------------------------------------------------------------
const P = 17;
const A = 0;
const B = 7;
const G = { x: 1, y: 5 }; // Base Point (Generator of order 9)

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
    if (P1.x === P2.x && mod(P1.y + P2.y, P) === 0) return null; // Point at infinity

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
    // 'E' is index 4
    let guessedKey = mod(maxIdx - 4, 26);
    return guessedKey;
}

// ----------------------------------------------------------------------------
// 4. INTERACTIVE CONTROLLER & UI BINDINGS
// ----------------------------------------------------------------------------
document.addEventListener('DOMContentLoaded', () => {
    // --- 4.1 Toy ECC Scheme Setup ---
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
        // Alice computes S = k * Q
        // Bob computes S = d * C1
        const S_Alice = scalarMult(kAlice, Q);
        const S_Bob = scalarMult(dBob, C1);
        sharedSecretDisplay.textContent = `${formatPoint(S_Alice)} (Shared)`;

        // 4. Hybrid Encryption: Mask characters with derived scalar from S
        // Key mask K = ((x_S * 7 + y_S * 11) mod 26) + 1
        const mask = mod(S_Alice.x * 7 + S_Alice.y * 11, 25) + 1;
        const ciphertext = caesarEncrypt(plaintext, mask);
        eccCipherDisplay.textContent = `${ciphertext}  [C1 = ${formatPoint(C1)}, StreamMask = +${mask}]`;

        // 5. Bob's Decryption:
        // Bob derives mask from S_Bob
        const bobMask = mod(S_Bob.x * 7 + S_Bob.y * 11, 25) + 1;
        const decrypted = caesarDecrypt(ciphertext, bobMask);
        eccDecryptedDisplay.textContent = decrypted;

        // Trace log
        eccTraceLog.innerHTML = `
[STEP 1: KEYGEN]
- Generator Point G = (${G.x}, ${G.y}) over F_17
- Bob Private Key d_B = ${dBob}
- Bob Public Key Q_B = d_B * G = ${formatPoint(Q)}

[STEP 2: ENCAPSULATION by Alice]
- Alice chooses ephemeral nonce k = ${kAlice}
- Alice Ephemeral Public Point C_1 = k * G = ${formatPoint(C1)}
- Alice Shared Secret Point S_A = k * Q_B = ${formatPoint(S_Alice)}
- Derived Symmetrical Mask = +${mask}

[STEP 3: ENCRYPTION]
- Plaintext: "${plaintext}"
- Ciphertext: "${ciphertext}"

[STEP 4: DECAPSULATION & DECRYPTION by Bob]
- Bob receives C_1 and computes S_B = d_B * C_1 = ${dBob} * ${formatPoint(C1)} = ${formatPoint(S_Bob)}
- Notice: S_B matches S_A identically due to the commutative scalar product!
- Bob unmasks ciphertext using -${bobMask} => "${decrypted}"
        `.trim();

        // Also synchronize with Caesar Cipher comparison
        document.getElementById('caesarPlaintext').value = plaintext;
        document.getElementById('caesarKey').value = mask;
        executeCaesar();
    }

    runEccBtn.addEventListener('click', executeToyECC);

    // --- 4.2 Caesar Cipher Section ---
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

    // --- 4.3 Statistical Cryptanalysis of Caesar Cipher ---
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
        const padding = 30;
        const chartWidth = width - 2 * padding;
        const chartHeight = height - 2 * padding;
        const barWidth = chartWidth / 26;

        // Find max percentage for scaling
        let maxPct = 0.16; // Standard max roughly 16% for 'E'
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

        // Draw Bars
        for (let i = 0; i < 26; i++) {
            const x = padding + i * barWidth;
            const letter = String.fromCharCode(65 + i);

            // Expected bar (gray outline / shaded)
            const expHeight = (ENGLISH_FREQ[i] / maxPct) * chartHeight;
            ctx.fillStyle = 'rgba(52, 152, 219, 0.25)';
            ctx.fillRect(x + 2, height - padding - expHeight, barWidth - 4, expHeight);
            ctx.strokeStyle = '#2980b9';
            ctx.strokeRect(x + 2, height - padding - expHeight, barWidth - 4, expHeight);

            // Observed bar (red/orange)
            const obsPct = totalLetters > 0 ? (observedCounts[i] / totalLetters) : 0;
            const obsHeight = (obsPct / maxPct) * chartHeight;
            ctx.fillStyle = 'rgba(231, 76, 60, 0.7)';
            ctx.fillRect(x + barWidth * 0.2, height - padding - obsHeight, barWidth * 0.6, obsHeight);

            // Label letter
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
        ctx.fillText('Expected English Frequency', width - 220, 25);

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

        // Compute Chi-Square for all 26 candidate keys
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

        // Sort by Chi-Square ascending (lower is better fit)
        results.sort((a, b) => a.chiSquare - b.chiSquare);

        const best = results[0];
        const isCorrect = (best.shift === actualShift);

        bestKeyBadge.innerHTML = `
            Recovered Key: <strong>k = ${best.shift}</strong> 
            (Actual: k = ${actualShift}) &nbsp;|&nbsp; 
            Status: <span style="color:${isCorrect ? '#27ae60' : '#e74c3c'}; font-weight:bold;">
                ${isCorrect ? '✓ SUCCESSFUL RECOVERY' : '✗ FAILED (Small Sample Variance)'}
            </span> &nbsp;|&nbsp; Minimum &Chi;&sup2; = ${best.chiSquare.toFixed(2)}
        `;

        // Render Table (Top 5 candidate keys)
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

        // Draw frequency distribution chart for the winning candidate
        drawDistributionChart(best.counts, best.totalLetters, best.shift);
    }

    runAnalysisBtn.addEventListener('click', runAnalysis);
    cryptoShiftSelect.addEventListener('change', runAnalysis);

    // --- 4.4 Multi-Method Benchmark Comparison ---
    const runBenchmarkBtn = document.getElementById('runBenchmarkBtn');
    const benchmarkResultCard = document.getElementById('benchmarkResultCard');

    runBenchmarkBtn.addEventListener('click', () => {
        const rawText = cryptoTextInput.value;
        const actualShift = parseInt(cryptoShiftSelect.value, 10);
        const ciphertext = caesarEncrypt(rawText, actualShift);
        const cleanLen = ciphertext.replace(/[^A-Z]/g, '').length;

        // 1. Chi-Square
        let bestChiK = -1, minChi = Infinity;
        for (let k = 0; k < 26; k++) {
            let res = computeChiSquare(ciphertext, k);
            if (res.chiSquare < minChi) {
                minChi = res.chiSquare;
                bestChiK = k;
            }
        }

        // 2. Single Letter ('E' Heuristic)
        const bestSingleK = computeSingleLetterKey(ciphertext);

        // 3. Bigram Score
        let bestBigramK = -1, maxBigram = -Infinity;
        for (let k = 0; k < 26; k++) {
            let s = computeBigramScore(ciphertext, k);
            if (s > maxBigram) {
                maxBigram = s;
                bestBigramK = k;
            }
        }

        benchmarkResultCard.innerHTML = `
            <div style="display:grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin-top: 1rem;">
                <div style="background:#fff; padding:1rem; border-radius:4px; border-left:4px solid ${bestChiK === actualShift ? '#27ae60' : '#e74c3c'};">
                    <strong>1. &Chi;&sup2; Frequency Vector:</strong><br>
                    Guessed Key: <strong>k = ${bestChiK}</strong><br>
                    Result: ${bestChiK === actualShift ? '<span style="color:#27ae60;">SUCCESS</span>' : '<span style="color:#e74c3c;">FAILED</span>'}
                </div>
                <div style="background:#fff; padding:1rem; border-radius:4px; border-left:4px solid ${bestSingleK === actualShift ? '#27ae60' : '#e74c3c'};">
                    <strong>2. Single-Letter ('E') Heuristic:</strong><br>
                    Guessed Key: <strong>k = ${bestSingleK}</strong><br>
                    Result: ${bestSingleK === actualShift ? '<span style="color:#27ae60;">SUCCESS</span>' : '<span style="color:#e74c3c;">FAILED</span>'}
                </div>
                <div style="background:#fff; padding:1rem; border-radius:4px; border-left:4px solid ${bestBigramK === actualShift ? '#27ae60' : '#e74c3c'};">
                    <strong>3. Bigram Log-Likelihood:</strong><br>
                    Guessed Key: <strong>k = ${bestBigramK}</strong><br>
                    Result: ${bestBigramK === actualShift ? '<span style="color:#27ae60;">SUCCESS</span>' : '<span style="color:#e74c3c;">FAILED</span>'}
                </div>
            </div>
            <p style="font-size:0.9rem; margin-top:0.8rem; color:#555;">
                <em>Statistical Insight:</em> For short texts ($N < 20$), the Single-Letter heuristic frequently fails because 'E' may not appear at all. Bigrams leverage pairwise Markov transitions, often breaking short ciphers where unigram &Chi;&sup2; exhibits ambiguous ties.
            </p>
        `;
    });

    // --- 4.5 Attack on Restricted-Space Toy ECC ---
    const runEccAttackBtn = document.getElementById('runEccAttackBtn');
    const eccAttackLog = document.getElementById('eccAttackLog');
    const eccAttackSummary = document.getElementById('eccAttackSummary');

    runEccAttackBtn.addEventListener('click', () => {
        const d_target = parseInt(bobPrivInput.value, 10);
        const Q_target = scalarMult(d_target, G);
        const startTime = performance.now();

        let iterations = 0;
        let recovered_d = null;
        let logLines = [];

        logLines.push(`[ATTACK INITIATED] Target Public Key Q = ${formatPoint(Q_target)}`);
        logLines.push(`[KEYSPACE EXPLOITATION] Restricted subgroup order n = 9 on F_17`);
        logLines.push(`Searching scalar space d in [1 .. 8]...`);

        for (let candidate = 1; candidate <= 8; candidate++) {
            iterations++;
            let testPoint = scalarMult(candidate, G);
            logLines.push(`Iteration #${candidate}: Testing scalar ${candidate} * G => ${formatPoint(testPoint)}`);

            if (testPoint && testPoint.x === Q_target.x && testPoint.y === Q_target.y) {
                recovered_d = candidate;
                logLines.push(`[MATCH FOUND!] Scalar multiplication matched target Public Key at d = ${recovered_d}!`);
                break;
            }
        }

        const elapsed = (performance.now() - startTime).toFixed(4);

        eccAttackSummary.innerHTML = `
            <div style="background:#e8f8f5; border:1px solid #27ae60; padding:1rem; border-radius:4px; margin-bottom:1rem;">
                <h4 style="margin:0 0 5px 0; color:#27ae60;">Private Key Recovered Successfully!</h4>
                Target Key: <strong>d = ${d_target}</strong> &nbsp;|&nbsp; 
                Recovered Key: <strong>d&#770; = ${recovered_d}</strong><br>
                Total Iterations: <strong>${iterations} / 8</strong> &nbsp;|&nbsp; 
                Time Elapsed: <strong>${elapsed} ms</strong>
            </div>
        `;

        eccAttackLog.textContent = logLines.join('\n');
    });

    // Initial triggers
    executeToyECC();
    runAnalysis();
});

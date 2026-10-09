// Parametri della Toy Curve: y^2 = x^3 + 7 (mod 17)
const P = 17;
const A = 0;
const B = 7;
const G = { x: 1, y: 5 }; // Punto Base (Generatore)

// Funzione Modulo sicura per JS (gestisce i numeri negativi)
function mod(n, p) {
    return ((n % p) + p) % p;
}

// Inverso Moltiplicativo Modulare (brute-force per p piccolo, andrebbe usato l'Euclideo Esteso)
function modInverse(a, m) {
    a = mod(a, m);
    for (let x = 1; x < m; x++) {
        if (mod(a * x, m) === 1) return x;
    }
    return null; // Non esiste
}

// Addizione di due punti sulla curva ellittica
function pointAdd(P1, P2) {
    // Gestione del Punto all'Infinito (rappresentato come null)
    if (!P1) return P2;
    if (!P2) return P1;

    // Se P1 == -P2, il risultato è il punto all'infinito
    if (P1.x === P2.x && mod(P1.y + P2.y, P) === 0) {
        return null; 
    }

    let m; // Pendenza (slope)
    if (P1.x === P2.x && P1.y === P2.y) {
        // Point Doubling (P1 + P1)
        // m = (3*x^2 + a) / (2*y) mod p
        let num = mod(3 * P1.x * P1.x + A, P);
        let den = modInverse(2 * P1.y, P);
        if (den === null) return null; // Dovrebbe accadere solo se y=0
        m = mod(num * den, P);
    } else {
        // Point Addition (P1 + P2)
        // m = (y2 - y1) / (x2 - x1) mod p
        let num = mod(P2.y - P1.y, P);
        let den = modInverse(P2.x - P1.x, P);
        if (den === null) return null;
        m = mod(num * den, P);
    }

    // x3 = m^2 - x1 - x2 mod p
    let x3 = mod(m * m - P1.x - P2.x, P);
    // y3 = m(x1 - x3) - y1 mod p
    let y3 = mod(m * (P1.x - x3) - P1.y, P);

    return { x: x3, y: y3 };
}

// Genera tutti i punti validi sulla curva per disegnarli in background
function getAllValidPoints() {
    let points = [];
    for (let x = 0; x < P; x++) {
        let rhs = mod(x * x * x + A * x + B, P);
        for (let y = 0; y < P; y++) {
            let lhs = mod(y * y, P);
            if (lhs === rhs) {
                points.push({ x, y });
            }
        }
    }
    return points;
}

// Visualizzazione "Guru Level" con Canvas
function drawVisualization(sequence) {
    const canvas = document.getElementById('curveCanvas');
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;
    const cellSize = width / P;

    ctx.clearRect(0, 0, width, height);

    // Disegna la griglia
    ctx.strokeStyle = '#ecf0f1';
    ctx.lineWidth = 1;
    for (let i = 0; i <= P; i++) {
        ctx.beginPath(); ctx.moveTo(i * cellSize, 0); ctx.lineTo(i * cellSize, height); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(0, i * cellSize); ctx.lineTo(width, i * cellSize); ctx.stroke();
    }

    // Disegna etichette assi (solo pari per pulizia)
    ctx.fillStyle = '#95a5a6';
    ctx.font = '10px sans-serif';
    for (let i = 0; i < P; i+=2) {
        ctx.fillText(i, i * cellSize + 2, height - 2); // Asse X
        ctx.fillText(P-1-i, 2, i * cellSize + 10);    // Asse Y (invertito per canvas)
    }

    // Coordinate canvas (0,0 in basso a sinistra)
    const toCanvas = (pt) => ({
        cx: pt.x * cellSize + cellSize / 2,
        cy: height - (pt.y * cellSize + cellSize / 2)
    });

    // 1. Disegna TUTTI i punti validi della curva (grigio chiaro)
    const allPoints = getAllValidPoints();
    ctx.fillStyle = '#bdc3c7';
    allPoints.forEach(pt => {
        let { cx, cy } = toCanvas(pt);
        ctx.beginPath(); ctx.arc(cx, cy, 4, 0, Math.PI * 2); ctx.fill();
    });

    if (!sequence || sequence.length === 0) return;

    // 2. Disegna le linee di connessione per evidenziare i salti G -> 2G -> 3G
    ctx.strokeStyle = 'rgba(52, 152, 219, 0.4)'; // Blu semitrasparente
    ctx.lineWidth = 2;
    ctx.beginPath();
    let start = toCanvas(sequence[0]);
    ctx.moveTo(start.cx, start.cy);
    for (let i = 1; i < sequence.length; i++) {
        if (!sequence[i]) break; // Punto all'infinito
        let next = toCanvas(sequence[i]);
        ctx.lineTo(next.cx, next.cy);
    }
    ctx.stroke();

    // 3. Disegna i punti della sequenza (blu scuro)
    ctx.fillStyle = '#2980b9';
    sequence.forEach((pt, index) => {
        if(!pt) return;
        let { cx, cy } = toCanvas(pt);
        ctx.beginPath(); ctx.arc(cx, cy, 6, 0, Math.PI * 2); ctx.fill();
        
        // Testo del multiplo (es. 2G)
        ctx.fillStyle = '#2c3e50';
        ctx.font = '12px bold sans-serif';
        ctx.fillText(`${index+1}G`, cx + 8, cy - 8);
        ctx.fillStyle = '#2980b9'; // reset per il prossimo punto
    });

    // 4. Evidenzia la Public Key Finale (rosso)
    let finalPoint = sequence[sequence.length - 1];
    if (finalPoint) {
        let { cx, cy } = toCanvas(finalPoint);
        ctx.fillStyle = '#e74c3c';
        ctx.beginPath(); ctx.arc(cx, cy, 8, 0, Math.PI * 2); ctx.fill();
    }
}

// Esecuzione al click
document.getElementById('computeBtn').addEventListener('click', () => {
    let k = parseInt(document.getElementById('privateKeyInput').value);
    if (isNaN(k) || k < 1) {
        alert("Please enter a positive integer for k.");
        return;
    }

    let logHtml = "";
    let sequence = [];
    let currentPoint = G;
    
    logHtml += `Starting computation on curve y² = x³ + 7 (mod 17)<br>`;
    logHtml += `G = (1, 5)<br><br>`;
    
    sequence.push(currentPoint);
    logHtml += `1G = (${currentPoint.x}, ${currentPoint.y})<br>`;

    for (let i = 2; i <= k; i++) {
        currentPoint = pointAdd(currentPoint, G);
        sequence.push(currentPoint);
        
        if (currentPoint === null) {
            logHtml += `${i}G = Point at Infinity (O)<br>`;
            break; // We reached the end of the group cycle
        } else {
            logHtml += `${i}G = (${currentPoint.x}, ${currentPoint.y})<br>`;
        }
    }

    if (currentPoint !== null) {
        logHtml += `<br><span class="highlight-key">PUBLIC KEY (Q) = (${currentPoint.x}, ${currentPoint.y})</span>`;
    }

    document.getElementById('logBox').innerHTML = logHtml;
    drawVisualization(sequence);
});

// Disegna il background al caricamento della pagina
drawVisualization([]);

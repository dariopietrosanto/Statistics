import { keccak_256 } from 'https://esm.sh/@noble/hashes@1.3.3/sha3';

function bytesToHex(bytes) {
    return Array.from(bytes).map(b => b.toString(16).padStart(2, '0')).join('');
}

function stringToBytes(str) {
    return new TextEncoder().encode(str);
}

document.getElementById('checksumBtn')?.addEventListener('click', () => {
    let rawInput = document.getElementById('rawAddress').value.trim().toLowerCase();
    
    // Rimuovi 0x se presente (non dovrebbe esserci visto il layout ma per sicurezza)
    if (rawInput.startsWith('0x')) {
        rawInput = rawInput.substring(2);
    }
    
    // Controlli base
    if (!/^[0-9a-f]{40}$/i.test(rawInput)) {
        alert("Please enter a valid 40-character hex address.");
        return;
    }
    
    // 1. Convert to lowercase
    const addressLower = rawInput.toLowerCase();
    
    // 2. Hash the lowercase address string using Keccak-256
    const hashBytes = keccak_256(stringToBytes(addressLower));
    const hashHex = bytesToHex(hashBytes);
    
    // 3. Apply EIP-55
    let checksummed = "";
    for (let i = 0; i < addressLower.length; i++) {
        const char = addressLower[i];
        // Se è un numero, rimane numero
        if (/[0-9]/.test(char)) {
            checksummed += char;
        } else {
            // Se è una lettera, controlliamo il nibble corrispondente nell'hash
            // parseInt(char, 16) se l'n-esimo nibble dell'hash è >= 8 -> MAIUSCOLA
            const hashNibble = parseInt(hashHex[i], 16);
            if (hashNibble >= 8) {
                checksummed += char.toUpperCase();
            } else {
                checksummed += char.toLowerCase();
            }
        }
    }
    
    // Update UI
    document.getElementById('demoOutput').style.display = 'block';
    document.getElementById('outLower').textContent = "0x" + addressLower;
    document.getElementById('outHash').textContent = hashHex;
    document.getElementById('outFinal').innerHTML = ""; // reset
    
    // Aggiungiamo formattazione visiva per le maiuscole modificate da EIP-55
    let highlightedFinal = "";
    for (let i = 0; i < checksummed.length; i++) {
        if (/[A-Z]/.test(checksummed[i])) {
            highlightedFinal += `<span style="color: #c0392b;">${checksummed[i]}</span>`; // evidenzia le maiuscole in rosso scuro
        } else {
            highlightedFinal += checksummed[i];
        }
    }
    document.getElementById('outFinal').innerHTML = highlightedFinal;
});

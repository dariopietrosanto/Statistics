const fs = require('fs');
const code = fs.readFileSync('hw1/script.js', 'utf8');

const domMock = `
const document = {
    getElementById: (id) => ({
        addEventListener: () => {},
        appendChild: () => {},
        value: '{"x":1,"y":5}',
        style: {},
        getContext: () => ({
            clearRect: () => {},
            beginPath: () => {},
            moveTo: () => {},
            lineTo: () => {},
            stroke: () => {},
            fillText: () => {},
            arc: () => {},
            fill: () => {}
        })
    }),
    createElement: () => ({})
};
`;

eval(domMock + code + "\nconsole.log('Math is OK! 2G = ', pointAdd(G, G));");

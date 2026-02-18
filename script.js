/* Interactive formula definitions and UI handling (expanded catalog)
     This file restores previously removed formulas (intersection solvers, sector, resistors, etc.)
     and adds a number of common bidirectional formulas.
*/

const formulaCatalog = {
    triangle: [
        {
            id: 'area',
            name: 'Area (A = base × height / 2)',
            params: [
                {key:'base', label:'Base (b)', type:'number'},
                {key:'height', label:'Height (h)', type:'number'},
                {key:'area', label:'Area (A)', type:'number'}
            ],
            solve: (vals)=>{
                const b = isFinite(vals.base)?vals.base:undefined;
                const h = isFinite(vals.height)?vals.height:undefined;
                const A = isFinite(vals.area)?vals.area:undefined;
                const known = [b,h,A].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two values (base, height, or area).');
                if(A===undefined) return {area: (b*h)/2};
                if(b===undefined) return {base: (2*A)/h};
                if(h===undefined) return {height: (2*A)/b};
                return {area:A, base:b, height:h};
            }
        },
        {
            id: 'pythagorean',
            name: 'Pythagorean (a, b, c)',
            params: [
                {key:'a', label:'Side a', type:'number'},
                {key:'b', label:'Side b', type:'number'},
                {key:'c', label:'Hypotenuse c', type:'number'}
            ],
            solve: (vals)=>{
                const a = isFinite(vals.a)?vals.a:undefined;
                const b = isFinite(vals.b)?vals.b:undefined;
                const c = isFinite(vals.c)?vals.c:undefined;
                const known = [a,b,c].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two values (a, b, or c).');
                if(c===undefined) return {c: Math.sqrt(a*a + b*b)};
                if(a===undefined) return {a: Math.sqrt(c*c - b*b)};
                if(b===undefined) return {b: Math.sqrt(c*c - a*a)};
                return {a,b,c};
            }
        },
        {
            id: 'sine rule', name: 'Sine Rule (sinA/a = sinB/b = sinC/c)',
            params: [
                {key:'a', label:'Side a', type:'number'},
                {key:'A', label:'Angle A (degrees)', type:'number'},
                {key:'b', label:'Side b', type:'number'},
                {key:'B', label:'Angle B (degrees)', type:'number'},
            ],
            solve: (vals)=>{
                const a = isFinite(vals.a)?vals.a:undefined;
                const A = isFinite(vals.A)?vals.A:undefined;
                const b = isFinite(vals.b)?vals.b:undefined;
                const B = isFinite(vals.B)?vals.B:undefined;
                const known = [a,A,b,B].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two values (a, A, b, or B).');
                const toRad = (deg) => deg * Math.PI / 180;
                if(a===undefined) return {a: b * Math.sin(toRad(A)) / Math.sin(toRad(B))};
                if(A===undefined) return {A: Math.asin(a * Math.sin(toRad(B)) / b) * 180 / Math.PI};
                if(b===undefined) return {b: a * Math.sin(toRad(B)) / Math.sin(toRad(A))};
                if(B===undefined) return {B: Math.asin(b * Math.sin(toRad(A)) / a) * 180 / Math.PI};
                return {a,A,b,B};
            }
        },
        {
            id: 'cosine rule', name: 'Cosine Rule (c² = a² + b² - 2ab cosC)',
            params: [
                {key:'a', label:'Side a', type:'number'},
                {key:'b', label:'Side b', type:'number'},
                {key:'c', label:'Side c', type:'number'},
                {key:'C', label:'Angle C (degrees)', type:'number'}
            ],
            solve: (vals)=>{
                const a = isFinite(vals.a)?vals.a:undefined;
                const b = isFinite(vals.b)?vals.b:undefined;
                const c = isFinite(vals.c)?vals.c:undefined;
                const C = isFinite(vals.C)?vals.C:undefined;
                const known = [a,b,c,C].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two values (a, b, c, or C).');
                const toRad = (deg) => deg * Math.PI / 180;
                if(c===undefined) return {c: Math.sqrt(a*a + b*b - 2*a*b*Math.cos(toRad(C)))};
                if(C===undefined) return {C: Math.acos((a*a + b*b - c*c) / (2*a*b)) * 180 / Math.PI};
                if(a===undefined) return {a: Math.sqrt(c*c + b*b - 2*b*c*Math.cos(toRad(C)))};
                if(b===undefined) return {b: Math.sqrt(c*c + a*a - 2*a*c*Math.cos(toRad(C)))};
                return {a,b,c,C};
            }
        }
    ],

    rectangle: [
        {
            id: 'area', name: 'Area (A = w × h)',
            params: [
                {key:'w', label:'Width (w)', type:'number'},
                {key:'h', label:'Height (h)', type:'number'},
                {key:'area', label:'Area (A)', type:'number'}
            ],
            solve: (vals)=>{
                const w = isFinite(vals.w)?vals.w:undefined;
                const h = isFinite(vals.h)?vals.h:undefined;
                const A = isFinite(vals.area)?vals.area:undefined;
                const known = [w,h,A].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two values (w, h, or area).');
                if(A===undefined) return {area: w*h};
                if(w===undefined) return {w: A/h};
                if(h===undefined) return {h: A/w};
                return {w,h,area:A};
            }
        },
        {
            id: 'perimeter', name: 'Perimeter (P = 2×(w+h))',
            params: [
                {key:'w', label:'Width (w)', type:'number'},
                {key:'h', label:'Height (h)', type:'number'},
                {key:'perimeter', label:'Perimeter (P)', type:'number'}
            ],
            solve: (vals)=>{
                const w = isFinite(vals.w)?vals.w:undefined;
                const h = isFinite(vals.h)?vals.h:undefined;
                const P = isFinite(vals.perimeter)?vals.perimeter:undefined;
                const known = [w,h,P].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two values (w, h, or perimeter).');
                if(P===undefined) return {perimeter: 2*(w+h)};
                if(w===undefined) return {w: P/2 - h};
                if(h===undefined) return {h: P/2 - w};
                return {w,h,perimeter:P};
            }
        }
    ],

    circle: [
        {
            id: 'area', name: 'Area (A = π r²)',
            params: [
                {key:'r', label:'Radius (r)', type:'number'},
                {key:'area', label:'Area (A)', type:'number'}
            ],
            solve: (vals)=>{
                const r = isFinite(vals.r)?vals.r:undefined;
                const A = isFinite(vals.area)?vals.area:undefined;
                const known = [r,A].filter(x=>x!==undefined).length;
                if(known < 1) throw new Error('Provide radius or area.');
                if(A===undefined) return {area: Math.PI * r * r};
                if(r===undefined) return {r: Math.sqrt(A / Math.PI)};
                return {r,A};
            }
        },
        {
            id: 'circumference', name: 'Circumference (C = 2π r)',
            params: [
                {key:'r', label:'Radius (r)', type:'number'},
                {key:'circumference', label:'Circumference (C)', type:'number'}
            ],
            solve: (vals)=>{
                const r = isFinite(vals.r)?vals.r:undefined;
                const C = isFinite(vals.circumference)?vals.circumference:undefined;
                const known = [r,C].filter(x=>x!==undefined).length;
                if(known < 1) throw new Error('Provide radius or circumference.');
                if(C===undefined) return {circumference: 2 * Math.PI * r};
                if(r===undefined) return {r: C / (2 * Math.PI)};
                return {r,C};
            }
        },
        {
            id: 'sector_area', name: 'Area of Sector (A = θ/360 × π r²)',
            params: [
                {key:'r', label:'Radius (r)', type:'number'},
                {key:'theta', label:'Central Angle θ (degrees)', type:'number'},
                {key:'area_sector', label:'Area of Sector', type:'number'}
            ],
            solve: (vals)=>{
                const r = isFinite(vals.r)?vals.r:undefined;
                const theta = isFinite(vals.theta)?vals.theta:undefined;
                const A_sector = isFinite(vals.area_sector)?vals.area_sector:undefined;
                const known = [r,theta,A_sector].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two values (r, θ, or area of sector).');
                if(A_sector===undefined) return {area_sector: (theta/360) * Math.PI * r * r};
                if(r===undefined) return {r: Math.sqrt((A_sector * 360) / (Math.PI * theta))};
                if(theta===undefined) return {theta: (A_sector * 360) / (Math.PI * r * r)};
                return {r,theta,A_sector};
            }
        }
    ],

    intersections: [
        {
            id: 'line_line', name: 'Line ↔ Line (a x + b y + c = 0)',
            params: [
                {key:'a1', label:'Line1 a', type:'number'},
                {key:'b1', label:'Line1 b', type:'number'},
                {key:'c1', label:'Line1 c', type:'number'},
                {key:'a2', label:'Line2 a', type:'number'},
                {key:'b2', label:'Line2 b', type:'number'},
                {key:'c2', label:'Line2 c', type:'number'}
            ],
            solve: (vals)=>{
                const a1 = isFinite(vals.a1)?vals.a1:0;
                const b1 = isFinite(vals.b1)?vals.b1:0;
                const c1 = isFinite(vals.c1)?vals.c1:0;
                const a2 = isFinite(vals.a2)?vals.a2:0;
                const b2 = isFinite(vals.b2)?vals.b2:0;
                const c2 = isFinite(vals.c2)?vals.c2:0;
                const D = a1*b2 - a2*b1;
                if(Math.abs(D) < 1e-12) throw new Error('Lines are parallel or coincident (no single intersection).');
                const x = (c2*b1 - c1*b2) / D;
                const y = (a2*c1 - a1*c2) / D;
                return {x: x, y: y};
            }
        },
        {
            id: 'line_parabola', name: 'Line (y = m x + b) ↔ Parabola (y = a x² + b x + c)',
            params: [
                {key:'p_a', label:'Parabola a', type:'number'},
                {key:'p_b', label:'Parabola b', type:'number'},
                {key:'p_c', label:'Parabola c', type:'number'},
                {key:'l_m', label:'Line slope m', type:'number'},
                {key:'l_b', label:'Line intercept b', type:'number'}
            ],
            solve: (vals)=>{
                const A = isFinite(vals.p_a)?vals.p_a:0;
                const B = isFinite(vals.p_b)?vals.p_b:0;
                const C = isFinite(vals.p_c)?vals.p_c:0;
                const m = isFinite(vals.l_m)?vals.l_m:0;
                const b = isFinite(vals.l_b)?vals.l_b:0;
                const a = A;
                const bb = B - m;
                const cc = C - b;
                if(Math.abs(a) < 1e-15){
                    if(Math.abs(bb) < 1e-15) throw new Error('No intersection or infinite intersections (degenerate).');
                    const x = -cc / bb;
                    const y = m * x + b;
                    return {x1: x, y1: y};
                }
                const disc = bb*bb - 4*a*cc;
                if(disc < 0) throw new Error('No real intersections (discriminant < 0).');
                const sqrtD = Math.sqrt(disc);
                const x1 = (-bb + sqrtD) / (2*a);
                const x2 = (-bb - sqrtD) / (2*a);
                const y1 = m * x1 + b;
                const y2 = m * x2 + b;
                if(Math.abs(x1 - x2) < 1e-12) return {x1: x1, y1: y1};
                return {x1: x1, y1: y1, x2: x2, y2: y2};
            }
        },
        {
            id: 'line_circle', name: 'Line ↔ Circle',
            params: [
                {key:'c_h', label:'Circle center h', type:'number'},
                {key:'c_k', label:'Circle center k', type:'number'},
                {key:'c_r', label:'Circle radius r', type:'number'},
                {key:'l_m', label:'Line slope m', type:'number'},
                {key:'l_b', label:'Line intercept b', type:'number'}
            ],
            solve: (vals)=>{
                const h = isFinite(vals.c_h)?vals.c_h:0;
                const k = isFinite(vals.c_k)?vals.c_k:0;
                const r = isFinite(vals.c_r)?vals.c_r:0;
                const m = isFinite(vals.l_m)?vals.l_m:0;
                const b = isFinite(vals.l_b)?vals.l_b:0;
                const a = 1 + m*m;
                const bb = -2*h + 2*m*(b - k);
                const cc = h*h + (b - k)*(b - k) - r*r;
                const disc = bb*bb - 4*a*cc;
                if(disc < 0) throw new Error('No real intersections (discriminant < 0).');
                const sqrtD = Math.sqrt(disc);
                const x1 = (-bb + sqrtD) / (2*a);
                const x2 = (-bb - sqrtD) / (2*a);
                const y1 = m*x1 + b;
                const y2 = m*x2 + b;
                if(Math.abs(x1 - x2) < 1e-12) return {x1:x1,y1:y1};
                return {x1:x1,y1:y1,x2:x2,y2:y2};
            }
        },
        {
            id: 'parabola_parabola', name: 'Parabola ↔ Parabola',
            params: [
                {key:'a1', label:'Parabola1 a', type:'number'},
                {key:'b1', label:'Parabola1 b', type:'number'},
                {key:'c1', label:'Parabola1 c', type:'number'},
                {key:'a2', label:'Parabola2 a', type:'number'},
                {key:'b2', label:'Parabola2 b', type:'number'},
                {key:'c2', label:'Parabola2 c', type:'number'}
            ],
            solve: (vals)=>{
                const a1 = isFinite(vals.a1)?vals.a1:0;
                const b1 = isFinite(vals.b1)?vals.b1:0;
                const c1 = isFinite(vals.c1)?vals.c1:0;
                const a2 = isFinite(vals.a2)?vals.a2:0;
                const b2 = isFinite(vals.b2)?vals.b2:0;
                const c2 = isFinite(vals.c2)?vals.c2:0;
                const a = a1 - a2;
                const bb = b1 - b2;
                const cc = c1 - c2;
                if(Math.abs(a) < 1e-15){
                    if(Math.abs(bb) < 1e-15) throw new Error('Parabolas are identical or parallel (no finite intersections).');
                    const x = -cc / bb;
                    const y = a1*x*x + b1*x + c1;
                    return {x1:x,y1:y};
                }
                const disc = bb*bb - 4*a*cc;
                if(disc < 0) throw new Error('No real intersections (discriminant < 0).');
                const sqrtD = Math.sqrt(disc);
                const x1 = (-bb + sqrtD) / (2*a);
                const x2 = (-bb - sqrtD) / (2*a);
                const y1 = a1*x1*x1 + b1*x1 + c1;
                const y2 = a1*x2*x2 + b1*x2 + c1;
                if(Math.abs(x1 - x2) < 1e-12) return {x1:x1,y1:y1};
                return {x1:x1,y1:y1,x2:x2,y2:y2};
            }
        },
        {
            id: 'circle_circle', name: 'Circle ↔ Circle',
            params: [
                {key:'x1', label:'Circle1 center x', type:'number'},
                {key:'y1', label:'Circle1 center y', type:'number'},
                {key:'r1', label:'Circle1 radius r', type:'number'},
                {key:'x2', label:'Circle2 center x', type:'number'},
                {key:'y2', label:'Circle2 center y', type:'number'},
                {key:'r2', label:'Circle2 radius r', type:'number'}
            ],
            solve: (vals)=>{
                const x1 = isFinite(vals.x1)?vals.x1:0;
                const y1 = isFinite(vals.y1)?vals.y1:0;
                const r1 = isFinite(vals.r1)?vals.r1:0;
                const x2 = isFinite(vals.x2)?vals.x2:0;
                const y2 = isFinite(vals.y2)?vals.y2:0;
                const r2 = isFinite(vals.r2)?vals.r2:0;
                const dx = x2 - x1;
                const dy = y2 - y1;
                const d = Math.hypot(dx, dy);
                if(d < 1e-15 && Math.abs(r1 - r2) < 1e-12) throw new Error('Circles are coincident (infinite intersections).');
                if(d > r1 + r2) throw new Error('Circles are separate (no intersections).');
                if(d < Math.abs(r1 - r2)) throw new Error('One circle contained within the other (no intersections).');
                const a = (r1*r1 - r2*r2 + d*d) / (2*d);
                const h2 = r1*r1 - a*a;
                const h = h2 > 0 ? Math.sqrt(h2) : 0;
                const xm = x1 + a * (dx) / d;
                const ym = y1 + a * (dy) / d;
                if(h === 0){
                    return {x1: xm, y1: ym};
                }
                const rx = -dy * (h / d);
                const ry = dx * (h / d);
                const xi1 = xm + rx;
                const yi1 = ym + ry;
                const xi2 = xm - rx;
                const yi2 = ym - ry;
                return {x1: xi1, y1: yi1, x2: xi2, y2: yi2};
            }
        }
    ],

    electrical: [
        {
            id: 'ohms_law', name: "Ohm's Law (V = I × R)",
            params: [
                {key:'v', label:'Voltage (V)', type:'number'},
                {key:'i', label:'Current (I)', type:'number'},
                {key:'r', label:'Resistance (R)', type:'number'}
            ],
            solve: (vals)=>{
                const V = isFinite(vals.v)?vals.v:undefined;
                const I = isFinite(vals.i)?vals.i:undefined;
                const R = isFinite(vals.r)?vals.r:undefined;
                const known = [V,I,R].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of V, I, R to compute the third.');
                if(R===undefined) return {r: V / I};
                if(I===undefined) return {i: V / R};
                if(V===undefined) return {v: I * R};
                return {v:V,i:I,r:R};
            }
        },
        {
            id: 'power', name: 'Power (P = V × I)',
            params: [
                {key:'p', label:'Power (P)', type:'number'},
                {key:'v', label:'Voltage (V)', type:'number'},
                {key:'i', label:'Current (I)', type:'number'}
            ],
            solve: (vals)=>{
                const P = isFinite(vals.p)?vals.p:undefined;
                const V = isFinite(vals.v)?vals.v:undefined;
                const I = isFinite(vals.i)?vals.i:undefined;
                const known = [P,V,I].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of P, V, I to compute the third.');
                if(P===undefined) return {p: V * I};
                if(V===undefined) return {v: P / I};
                if(I===undefined) return {i: P / V};
                return {p:P,v:V,i:I};
            }
        },
        {
            id: 'res_series', name: 'Resistors (series) R = R1 + R2',
            params: [
                {key:'R', label:'Total R', type:'number'},
                {key:'R1', label:'R1', type:'number'},
                {key:'R2', label:'R2', type:'number'}
            ],
            solve: (vals)=>{
                const R = isFinite(vals.R)?vals.R:undefined;
                const R1 = isFinite(vals.R1)?vals.R1:undefined;
                const R2 = isFinite(vals.R2)?vals.R2:undefined;
                const known = [R,R1,R2].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of R, R1, R2.');
                if(R===undefined) return {R: R1 + R2};
                if(R1===undefined) return {R1: R - R2};
                if(R2===undefined) return {R2: R - R1};
                return {R,R1,R2};
            }
        },
        {
            id: 'res_parallel', name: 'Resistors (parallel) 1/R = 1/R1 + 1/R2',
            params: [
                {key:'R', label:'Total R', type:'number'},
                {key:'R1', label:'R1', type:'number'},
                {key:'R2', label:'R2', type:'number'}
            ],
            solve: (vals)=>{
                const R = isFinite(vals.R)?vals.R:undefined;
                const R1 = isFinite(vals.R1)?vals.R1:undefined;
                const R2 = isFinite(vals.R2)?vals.R2:undefined;
                const known = [R,R1,R2].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of R, R1, R2.');
                if(R===undefined) return {R: 1 / (1 / R1 + 1 / R2)};
                if(R1===undefined) return {R1: 1 / (1 / R - 1 / R2)};
                if(R2===undefined) return {R2: 1 / (1 / R - 1 / R1)};
                return {R,R1,R2};
            }
        }
    ],

    mechanics: [
        {
            id: 'force', name: 'Force (F = m × a)',
            params: [
                {key:'F', label:'Force (F)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'a', label:'Acceleration (a)', type:'number'}
            ],
            solve: (vals)=>{
                const F = isFinite(vals.F)?vals.F:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;
                const a = isFinite(vals.a)?vals.a:undefined;
                const known = [F,m,a].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of F, m, a.');
                if(F===undefined) return {F: m * a};
                if(m===undefined) return {m: F / a};
                if(a===undefined) return {a: F / m};
                return {F,m,a};
            }
        },
        {
            id: 'momentum', name: 'Momentum (p = m × v)',
            params: [
                {key:'p', label:'Momentum (p)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'v', label:'Velocity (v)', type:'number'}
            ],
            solve: (vals)=>{
                const p = isFinite(vals.p)?vals.p:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;
                const v = isFinite(vals.v)?vals.v:undefined;
                const known = [p,m,v].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of p, m, v.');
                if(p===undefined) return {p: m * v};
                if(m===undefined) return {m: p / v};
                if(v===undefined) return {v: p / m};
                return {p,m,v};
            }
        },
        {
            id: 'work', name: 'Work (W = F × d)',
            params: [
                {key:'W', label:'Work (W)', type:'number'},
                {key:'F', label:'Force (F)', type:'number'},
                {key:'d', label:'Distance (d)', type:'number'}
            ],
            solve: (vals)=>{
                const W = isFinite(vals.W)?vals.W:undefined;
                const F = isFinite(vals.F)?vals.F:undefined;
                const d = isFinite(vals.d)?vals.d:undefined;
                const known = [W,F,d].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of W, F, d.');
                if(W===undefined) return {W: F * d};
                if(F===undefined) return {F: W / d};
                if(d===undefined) return {d: W / F};
                return {W,F,d};
            }
        },
        {
            id: 'pressure', name: 'Pressure (P = F / A)',
            params: [
                {key:'P', label:'Pressure (P)', type:'number'},
                {key:'F', label:'Force (F)', type:'number'},
                {key:'A', label:'Area (A)', type:'number'}
            ],
            solve: (vals)=>{
                const P = isFinite(vals.P)?vals.P:undefined;
                const F = isFinite(vals.F)?vals.F:undefined;
                const A = isFinite(vals.A)?vals.A:undefined;
                const known = [P,F,A].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of P, F, A.');
                if(P===undefined) return {P: F / A};
                if(F===undefined) return {F: P * A};
                if(A===undefined) return {A: F / P};
                return {P,F,A};
            }
        },
        {
            id: 'density', name: 'Density (ρ = m / V)',
            params: [
                {key:'rho', label:'Density (ρ)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'V', label:'Volume (V)', type:'number'}
            ],
            solve: (vals)=>{
                const rho = isFinite(vals.rho)?vals.rho:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;    
                const V = isFinite(vals.V)?vals.V:undefined;
                const known = [rho,m,V].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of ρ, m, V.');
                if(rho===undefined) return {rho: m / V};
                if(m===undefined) return {m: rho * V};
                if(V===undefined) return {V: m / rho};
                return {rho,m,V};
            }
        },
        {
            id: 'force_direction', name: 'Force Direction (F = m × g)',
            params: [
                {key:'F', label:'Force (F)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'g', label:'Gravity (g) [default 9.81]', type:'number'}
            ],
            solve: (vals)=>{
                const F = isFinite(vals.F)?vals.F:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;
                const g = isFinite(vals.g)?vals.g:9.81;
                const known = [F,m].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two of F, m (g defaults to 9.81 if omitted).');
                if(F===undefined) return {F: m * g};
                if(m===undefined) return {m: F / g};
                return {F,m,g};
            }
        },
        {
            id: 'acceleration', name: 'Acceleration (a = Δv / Δt)',
            params: [
                {key:'a', label:'Acceleration (a)', type:'number'},
                {key:'delta_v', label:'Change in Velocity (Δv)', type:'number'},
                {key:'delta_t', label:'Change in Time (Δt)', type:'number'}
            ],
            solve: (vals)=>{
                const a = isFinite(vals.a)?vals.a:undefined;
                const delta_v = isFinite(vals.delta_v)?vals.delta_v:undefined;
                const delta_t = isFinite(vals.delta_t)?vals.delta_t:undefined;
                const known = [a,delta_v,delta_t].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of a, Δv, Δt.');
                if(a===undefined) return {a: delta_v / delta_t};
                if(delta_v===undefined) return {delta_v: a * delta_t};
                if(delta_t===undefined) return {delta_t: delta_v / a};
                return {a,delta_v,delta_t};
            }
        },
        {
            id: 'velocity', name: 'Velocity (v = d / t)',
            params: [
                {key:'v', label:'Velocity (v)', type:'number'},
                {key:'d', label:'Distance (d)', type:'number'},
                {key:'t', label:'Time (t)', type:'number'}
            ],
            solve: (vals)=>{
                const v = isFinite(vals.v)?vals.v:undefined;
                const d = isFinite(vals.d)?vals.d:undefined;
                const t = isFinite(vals.t)?vals.t:undefined;
                const known = [v,d,t].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of v, d, t.');
                if(v===undefined) return {v: d / t};
                if(d===undefined) return {d: v * t};
                if(t===undefined) return {t: d / v};
                return {v,d,t};
            }
        },
        {
            id: 'displacement', name: 'Displacement (s = ut + 1/2 at²)',
            params: [
                {key:'s', label:'Displacement (s)', type:'number'},
                {key:'u', label:'Initial Velocity (u)', type:'number'},
                {key:'a', label:'Acceleration (a)', type:'number'},
                {key:'t', label:'Time (t)', type:'number'}
            ],
            solve: (vals)=>{
                const s = isFinite(vals.s)?vals.s:undefined;
                const u = isFinite(vals.u)?vals.u:undefined;
                const a = isFinite(vals.a)?vals.a:undefined;
                const t = isFinite(vals.t)?vals.t:undefined;
                const known = [s,u,a,t].filter(x=>x!==undefined).length;
                if(known < 3) throw new Error('Provide at least three values out of s, u, a, t.');
                if(s===undefined) return {s: u * t + 0.5 * a * t * t};
                if(u===undefined) return {u: (s - 0.5 * a * t * t) / t};
                if(a===undefined) return {a: (2 * (s - u * t)) / (t * t)};
                if(t===undefined) {
                    const discriminant = u*u + 2*a*s;
                    if(discriminant < 0) throw new Error('No real solution for time (discriminant < 0).');
                    const sqrtD = Math.sqrt(discriminant);
                    const t1 = ( -u + sqrtD ) / a;
                    const t2 = ( -u - sqrtD ) / a;
                    return {t1: t1, t2: t2};
                }
                return {s,u,a,t};
            }
        },
        {
            id: 'hookes_law', name: "Hooke's Law (F = k × x)",
            params: [
                {key:'F', label:'Force (F)', type:'number'},
                {key:'k', label:'Spring Constant (k)', type:'number'},
                {key:'x', label:'Displacement (x)', type:'number'}
            ],
            solve: (vals)=>{
                const F = isFinite(vals.F)?vals.F:undefined;
                const k = isFinite(vals.k)?vals.k:undefined;
                const x = isFinite(vals.x)?vals.x:undefined;
                const known = [F,k,x].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of F, k, x.');
                if(F===undefined) return {F: k * x};
                if(k===undefined) return {k: F / x};
                if(x===undefined) return {x: F / k};
                return {F,k,x};
            }
        },
        {
            id: 'centripetal_force', name: 'Centripetal Force (F = m v² / r)',
            params: [
                {key:'F', label:'Centripetal Force (F)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'v', label:'Velocity (v)', type:'number'},
                {key:'r', label:'Radius (r)', type:'number'}
            ],
            solve: (vals)=>{
                const F = isFinite(vals.F)?vals.F:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;
                const v = isFinite(vals.v)?vals.v:undefined;
                const r = isFinite(vals.r)?vals.r:undefined;
                const known = [F,m,v,r].filter(x=>x!==undefined).length;
                if(known < 3) throw new Error('Provide at least three values out of F, m, v, r.');
                if(F===undefined) return {F: (m * v * v) / r};
                if(m===undefined) return {m: (F * r) / (v * v)};
                if(v===undefined) return {v: Math.sqrt((F * r) / m)};
                if(r===undefined) return {r: (m * v * v) / F};
                return {F,m,v,r};
            }
        }
    ],

    energy: [
        {
            id: 'kinetic', name: 'Kinetic Energy (KE = 1/2 m v²)',
            params: [
                {key:'ke', label:'Kinetic Energy (KE)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'v', label:'Velocity (v)', type:'number'}
            ],
            solve: (vals)=>{
                const ke = isFinite(vals.ke)?vals.ke:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;
                const v = isFinite(vals.v)?vals.v:undefined;
                const known = [ke,m,v].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of KE, m, v.');
                if(ke===undefined) return {ke: 0.5 * m * v * v};
                if(m===undefined) return {m: (2*ke) / (v*v)};
                if(v===undefined) return {v: Math.sqrt((2*ke)/m)};
                return {ke,m,v};
            }
        },
        {
            id: 'potential', name: 'Potential Energy (PE = m g h)',
            params: [
                {key:'pe', label:'Potential Energy (PE)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'g', label:'Gravity (g) [default 9.81]', type:'number'},
                {key:'h', label:'Height (h)', type:'number'}
            ],
            solve: (vals)=>{
                const pe = isFinite(vals.pe)?vals.pe:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;
                const g = isFinite(vals.g)?vals.g:9.81;
                const h = isFinite(vals.h)?vals.h:undefined;
                const known = [pe,m,h].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide at least two of PE, m, h (g defaults to 9.81 if omitted).');
                if(pe===undefined) return {pe: m * g * h};
                if(m===undefined) return {m: pe / (g * h)};
                if(h===undefined) return {h: pe / (m * g)};
                return {pe,m,g,h};
            }
        }
    ],

    waves: [
        {
            id: 'wave_speed', name: 'Wave Speed (v = f × λ)',
            params: [
                {key:'v', label:'Speed (v)', type:'number'},
                {key:'f', label:'Frequency (f)', type:'number'},
                {key:'lambda', label:'Wavelength (λ)', type:'number'}
            ],
            solve: (vals)=>{
                const v = isFinite(vals.v)?vals.v:undefined;
                const f = isFinite(vals.f)?vals.f:undefined;
                const lambda = isFinite(vals.lambda)?vals.lambda:undefined;
                const known = [v,f,lambda].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of v, f, λ.');
                if(v===undefined) return {v: f * lambda};
                if(f===undefined) return {f: v / lambda};
                if(lambda===undefined) return {lambda: v / f};
                return {v,f,lambda};
            }
        }
    ],

    optics: [
        {
            id: 'lens', name: 'Lens equation (1/f = 1/do + 1/di)',
            params: [
                {key:'f', label:'Focal length (f)', type:'number'},
                {key:'do', label:'Object distance (do)', type:'number'},
                {key:'di', label:'Image distance (di)', type:'number'}
            ],
            solve: (vals)=>{
                const f = isFinite(vals.f)?vals.f:undefined;
                const do_ = isFinite(vals.do)?vals.do:undefined;
                const di = isFinite(vals.di)?vals.di:undefined;
                const known = [f,do_,di].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of f, do, di.');
                if(f===undefined) return {f: 1 / (1/do_ + 1/di)};
                if(do_===undefined) return {do: 1 / (1/f - 1/di)};
                if(di===undefined) return {di: 1 / (1/f - 1/do_)};
                return {f,do:do_,di};
            }
        }
    ],

    finance: [
        {
            id: 'simple_interest', name: 'Simple Interest (A = P (1 + r t))',
            params: [
                {key:'A', label:'Amount (A)', type:'number'},
                {key:'P', label:'Principal (P)', type:'number'},
                {key:'r', label:'Rate (r, decimal)', type:'number'},
                {key:'t', label:'Time (t)', type:'number'}
            ],
            solve: (vals)=>{
                const A = isFinite(vals.A)?vals.A:undefined;
                const P = isFinite(vals.P)?vals.P:undefined;
                const r = isFinite(vals.r)?vals.r:undefined;
                const t = isFinite(vals.t)?vals.t:undefined;
                const known = [A,P,r,t].filter(x=>x!==undefined).length;
                if(known < 3) throw new Error('Provide three values out of A, P, r, t.');
                if(A===undefined) return {A: P * (1 + r * t)};
                if(P===undefined) return {P: A / (1 + r * t)};
                if(r===undefined) return {r: (A / P - 1) / t};
                if(t===undefined) return {t: (A / P - 1) / r};
                return {A,P,r,t};
            }
        }
    ],

    misc: [
        {
            id: 'density', name: 'Density (ρ = m / V)',
            params: [
                {key:'rho', label:'Density (ρ)', type:'number'},
                {key:'m', label:'Mass (m)', type:'number'},
                {key:'V', label:'Volume (V)', type:'number'}
            ],
            solve: (vals)=>{
                const rho = isFinite(vals.rho)?vals.rho:undefined;
                const m = isFinite(vals.m)?vals.m:undefined;
                const V = isFinite(vals.V)?vals.V:undefined;
                const known = [rho,m,V].filter(x=>x!==undefined).length;
                if(known < 2) throw new Error('Provide two values out of density, mass, volume.');
                if(rho===undefined) return {rho: m / V};
                if(m===undefined) return {m: rho * V};
                if(V===undefined) return {V: m / rho};
                return {rho,m,V};
            }
        },
        {
            id: 'frequency', name: 'Frequency / Period (f = 1 / T)',
            params: [
                {key:'f', label:'Frequency (f)', type:'number'},
                {key:'T', label:'Period (T)', type:'number'}
            ],
            solve: (vals)=>{
                const f = isFinite(vals.f)?vals.f:undefined;
                const T = isFinite(vals.T)?vals.T:undefined;
                if(f===undefined && T===undefined) throw new Error('Provide frequency or period.');
                if(f===undefined) return {f: 1 / T};
                if(T===undefined) return {T: 1 / f};
                return {f,T};
            }
        }
    ],
    converter: [
        {
            id: 'celsius_fahrenheit', name: 'Celsius ↔ Fahrenheit',
            params: [
                {key:'C', label:'Celsius (°C)', type:'number'},
                {key:'F', label:'Fahrenheit (°F)', type:'number'}
            ],
            solve: (vals)=>{
                const C = isFinite(vals.C)?vals.C:undefined;
                const F = isFinite(vals.F)?vals.F:undefined;
                if(C===undefined && F===undefined) throw new Error('Provide Celsius or Fahrenheit.');
                if(C===undefined) return {C: (F - 32) * 5 / 9};
                if(F===undefined) return {F: C * 9 / 5 + 32};
                return {C,F};
            }
        },
        {
            id: 'km_miles', name: 'Kilometers ↔ Miles',
            params: [
                {key:'km', label:'Kilometers (km)', type:'number'},
                {key:'miles', label:'Miles', type:'number'}
            ],
            solve: (vals)=>{
                const km = isFinite(vals.km)?vals.km:undefined;
                const miles = isFinite(vals.miles)?vals.miles:undefined;
                if(km===undefined && miles===undefined) throw new Error('Provide kilometers or miles.');
                if(km===undefined) return {km: miles * 1.60934};
                if(miles===undefined) return {miles: km / 1.60934};
                return {km,miles};
            }
        },
        {
            id: 'degrees_radians', name: 'Degrees ↔ Radians',
            params: [
                {key:'degrees', label:'Degrees (°)', type:'number'},
                {key:'radians', label:'Radians', type:'number'}
            ],
            solve: (vals)=>{
                const degrees = isFinite(vals.degrees)?vals.degrees:undefined;
                const radians = isFinite(vals.radians)?vals.radians:undefined;
                if(degrees===undefined && radians===undefined) throw new Error('Provide degrees or radians.');
                if(degrees===undefined) return {degrees: radians * (180 / Math.PI)};
                if(radians===undefined) return {radians: degrees * (Math.PI / 180)};
                return {degrees,radians};
            }
        },
        {
            id: 'kg_pounds', name: 'Kilograms ↔ Pounds',
            params: [
                {key:'kg', label:'Kilograms (kg)', type:'number'},  
                {key:'pounds', label:'Pounds (lbs)', type:'number'}
            ],
            solve: (vals)=>{
                const kg = isFinite(vals.kg)?vals.kg:undefined;
                const pounds = isFinite(vals.pounds)?vals.pounds:undefined;
                if(kg===undefined && pounds===undefined) throw new Error('Provide kilograms or pounds.');
                if(kg===undefined) return {kg: pounds * 0.453592};
                if(pounds===undefined) return {pounds: kg / 0.453592};
                return {kg,pounds};
            }
        },
        {
            id: 'liters_gallons', name: 'Liters ↔ Gallons',
            params: [
                {key:'liters', label:'Liters (L)', type:'number'},
                {key:'gallons', label:'Gallons (gal)', type:'number'}
            ],
            solve: (vals)=>{
                const liters = isFinite(vals.liters)?vals.liters:undefined;
                const gallons = isFinite(vals.gallons)?vals.gallons:undefined;
                if(liters===undefined && gallons===undefined) throw new Error('Provide liters or gallons.');
                if(liters===undefined) return {liters: gallons * 3.78541};
                if(gallons===undefined) return {gallons: liters / 3.78541};
                return {liters,gallons};
            }
        },
        {
            id: 'cubic_meters_cubic_feet', name: 'Cubic Meters ↔ Cubic Feet',
            params: [
                {key:'cubic_meters', label:'Cubic Meters (m³)', type:'number'},
                {key:'cubic_feet', label:'Cubic Feet (ft³)', type:'number'}
            ],
            solve: (vals)=>{
                const cubic_meters = isFinite(vals.cubic_meters)?vals.cubic_meters:undefined;
                const cubic_feet = isFinite(vals.cubic_feet)?vals.cubic_feet:undefined;
                if(cubic_meters===undefined && cubic_feet===undefined) throw new Error('Provide cubic meters or cubic feet.');
                if(cubic_meters===undefined) return {cubic_meters: cubic_feet * 0.0283168};
                if(cubic_feet===undefined) return {cubic_feet: cubic_meters / 0.0283168};
                return {cubic_meters,cubic_feet};
            }
        },
        {
            id: 'newtons_lbf', name: 'Newtons ↔ Pounds-force',
            params: [
                {key:'newtons', label:'Newtons (N)', type:'number'},
                {key:'lbf', label:'Pounds-force (lbf)', type:'number'}
            ],
            solve: (vals)=>{
                const newtons = isFinite(vals.newtons)?vals.newtons:undefined;
                const lbf = isFinite(vals.lbf)?vals.lbf:undefined;
                if(newtons===undefined && lbf===undefined) throw new Error('Provide newtons or pounds-force.');
                if(newtons===undefined) return {newtons: lbf * 4.44822};
                if(lbf===undefined) return {lbf: newtons / 4.44822};
                return {newtons,lbf};
            }
        },
        {
            id: 'watts_hp', name: 'Watts ↔ Horsepower',
            params: [
                {key:'watts', label:'Watts (W)', type:'number'},
                {key:'hp', label:'Horsepower (hp)', type:'number'}
            ],
            solve: (vals)=>{
                const watts = isFinite(vals.watts)?vals.watts:undefined;
                const hp = isFinite(vals.hp)?vals.hp:undefined;
                if(watts===undefined && hp===undefined) throw new Error('Provide watts or horsepower.');
                if(watts===undefined) return {watts: hp * 745.7};
                if(hp===undefined) return {hp: watts / 745.7};
                return {watts,hp};
            }
        },
        {
            id: 'joules_calories', name: 'Joules ↔ Calories',
            params: [
                {key:'joules', label:'Joules (J)', type:'number'},
                {key:'calories', label:'Calories (cal)', type:'number'}
            ],
            solve: (vals)=>{
                const joules = isFinite(vals.joules)?vals.joules:undefined;
                const calories = isFinite(vals.calories)?vals.calories:undefined;
                if(joules===undefined && calories===undefined) throw new Error('Provide joules or calories.');
                if(joules===undefined) return {joules: calories * 4.184};
                if(calories===undefined) return {calories: joules / 4.184};
                return {joules,calories};
            }
        }
    ]
};

// ----- Helpers and UI bindings -----
document.addEventListener('DOMContentLoaded', ()=>{
    const catSel = document.getElementById('formula_type');
    const formSel = document.getElementById('formula_name');
    const loadBtn = document.getElementById('load_formula');
    const clearBtn = document.getElementById('clear');
    const paramsContainer = document.getElementById('params');
    const answerDiv = document.getElementById('answer-display');

    // populate categories
    Object.keys(formulaCatalog).forEach(cat=>{
        const opt = document.createElement('option');
        opt.value = cat;
        opt.textContent = cat.charAt(0).toUpperCase() + cat.slice(1);
        catSel.appendChild(opt);
    });

    function populateFormulas(){
        formSel.innerHTML = '';
        const cat = catSel.value;
        if(!cat || !formulaCatalog[cat]) return;
        formulaCatalog[cat].forEach(f=>{
            const opt = document.createElement('option');
            opt.value = f.id;
            opt.textContent = f.name;
            formSel.appendChild(opt);
        });
    }

    // parse flexible numeric inputs: plain numbers, 1e3, or forms like 1.23*10^4
    function parseNumericInput(str){
        if(str === undefined || str === null) return undefined;
        const s = String(str).trim();
        if(s === '') return undefined;
        const direct = Number.parseFloat(s);
        if(!Number.isNaN(direct) && /^[+-]?[0-9]*\.?[0-9]+([eE][+-]?\d+)?$/.test(s)) return direct;
        const m = s.match(/^([+-]?[0-9]*\.?[0-9]+)\s*(?:\*|×)\s*10\s*(?:\^|\u02c6)\s*([+-]?\d+)$/);
        if(m){
            const mant = Number.parseFloat(m[1]);
            const exp = Number.parseInt(m[2],10);
            if(Number.isFinite(mant) && Number.isFinite(exp)) return mant * Math.pow(10, exp);
        }
        const m2 = s.match(/^([+-]?[0-9]*\.?[0-9]+)\s+10\s*(?:\^|\u02c6)\s*([+-]?\d+)$/);
        if(m2){
            const mant = Number.parseFloat(m2[1]);
            const exp = Number.parseInt(m2[2],10);
            if(Number.isFinite(mant) && Number.isFinite(exp)) return mant * Math.pow(10, exp);
        }
        return NaN;
    }

    function formatScientific(n){
        if(!Number.isFinite(n)) return String(n);
        if(n === 0) return {mantissa: 0, exponent: 0, str: '0 × 10^0'};
        const sign = n < 0 ? -1 : 1;
        const abs = Math.abs(n);
        const exponent = Math.floor(Math.log10(abs));
        const mantissa = sign * (abs / Math.pow(10, exponent));
        const mStr = Number.parseFloat(mantissa.toFixed(6)).toString();
        return {mantissa: mStr, exponent: exponent, str: `${mStr} × 10^${exponent}`};
    }

    // safe evaluator for expressions including common math functions and constants
    // Supported: + - * / ** ^ ( ) numbers, e/E, and functions: sin, cos, tan, asin, acos, atan, sqrt, ln, log, exp, abs
    // Constants: pi, e
    // Trig functions respect `degreeMode` (if true, inputs are degrees and outputs of inverse funcs are degrees)
    function safeEvalExpression(expr){
        if(expr === undefined || expr === null) return undefined;
        let s = String(expr).trim();
        if(s === '') return undefined;
        // normalize common symbols
        s = s.replace(/π/g, 'pi');
        s = s.replace(/\^/g, '**');

        // identify identifiers and ensure they are in the allowed whitelist
        // remove numeric literals first (including scientific notation) so 'e' inside 1e3 isn't considered an identifier
        const numericPattern = /(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?/g;
        const withoutNums = s.replace(numericPattern, ' ');
        const allowed = new Set(['sin','cos','tan','asin','acos','atan','sqrt','ln','log','exp','abs','pi','e','pow']);
        const ids = withoutNums.match(/[A-Za-z_][A-Za-z0-9_]*/g) || [];
        for(const id of ids){
            if(!allowed.has(id)) return NaN;
        }

        // build wrapper functions that respect degreeMode
        const toRad = (x) => (degreeMode ? x * Math.PI / 180 : x);
        const fromRad = (x) => (degreeMode ? x * 180 / Math.PI : x);
        const sin = (x) => Math.sin(toRad(x));
        const cos = (x) => Math.cos(toRad(x));
        const tan = (x) => Math.tan(toRad(x));
        const asin = (x) => fromRad(Math.asin(x));
        const acos = (x) => fromRad(Math.acos(x));
        const atan = (x) => fromRad(Math.atan(x));
        const sqrt = (x) => Math.sqrt(x);
        const ln = (x) => Math.log(x);
        const log = (x) => Math.log10 ? Math.log10(x) : Math.log(x) / Math.LN10;
        const exp = (x) => Math.exp(x);
        const abs = (x) => Math.abs(x);
        const pi = Math.PI;
        const e = Math.E;

        // allow only safe characters (digits, operators, parentheses, decimal, comma, identifiers handled above)
        if(!/^[0-9eE+\-*/().,\s*A-Za-z_\*\*]+$/.test(s)) return NaN;

        try{
            // construct Function with allowed names as parameters to avoid exposing globals
            const names = ['sin','cos','tan','asin','acos','atan','sqrt','ln','log','exp','abs','pi','e','pow'];
            const fn = new Function(...names, 'return (' + s + ')');
            const args = [sin, cos, tan, asin, acos, atan, sqrt, ln, log, exp, abs, pi, e, Math.pow];
            const res = fn(...args);
            if(typeof res === 'number' && Number.isFinite(res)) return res;
            return NaN;
        }catch(e){
            return NaN;
        }
    }

    catSel.addEventListener('change', ()=>{
        populateFormulas();
        paramsContainer.innerHTML = '';
        answerDiv.textContent = 'Result will appear here';
    });

    loadBtn.addEventListener('click', ()=>{
        const cat = catSel.value;
        const fid = formSel.value;
        if(!cat){ answerDiv.textContent = 'Select a category.'; return; }
        const formula = (formulaCatalog[cat] || []).find(x=>x.id===fid) || formulaCatalog[cat][0];
        if(!formula){ answerDiv.textContent = 'Select a formula.'; return; }

        paramsContainer.innerHTML = '';
        formula.params.forEach(p=>{
            const wrapper = document.createElement('div');
            wrapper.className = 'label_inbut';
            const label = document.createElement('label');
            label.textContent = p.label;
            label.htmlFor = `param_${p.key}`;
            const input = document.createElement('input');
            // use text input so users can type expressions (2+3, 1e3, etc.) and use on-screen calculator
            input.type = 'text';
            input.id = `param_${p.key}`;
            input.dataset.key = p.key;
            input.dataset.origType = p.type || 'number';
            input.value = '';
            input.placeholder = 'Leave blank to compute';
            wrapper.appendChild(label);
            wrapper.appendChild(input);
            paramsContainer.appendChild(wrapper);
        });

        const computeBtn = document.createElement('button');
        computeBtn.className = 'btn';
        computeBtn.textContent = 'Compute';
        // safeEvalExpression is declared in the outer scope so calculator can reuse it

        computeBtn.addEventListener('click', ()=>{
            const inputs = paramsContainer.querySelectorAll('input');
            const vals = {};
            let parseError = false;
            inputs.forEach(inp=>{
                const k = inp.dataset.key;
                const raw = inp.value.trim();
                // try direct numeric parse first
                let n = parseNumericInput(raw);
                if(Number.isNaN(n)){
                    // try evaluating an arithmetic expression
                    n = safeEvalExpression(raw);
                }
                if(Number.isNaN(n)) parseError = true;
                vals[k] = (n === undefined) ? undefined : n;
            });
            if(parseError){ answerDiv.textContent = 'One or more inputs could not be parsed as numbers or expressions.'; return; }
            try{
                const solved = formula.solve(vals);
                const allKeys = new Set([...Object.keys(vals), ...Object.keys(solved)]);
                const pairs = [];
                allKeys.forEach(k=>{
                    const rawVal = (k in solved && solved[k] !== undefined) ? solved[k] : vals[k];
                    if(rawVal === undefined){ pairs.push(`${k}: —`); return; }
                    if(typeof rawVal === 'number'){
                        const dec = Number.parseFloat(rawVal.toFixed(6));
                        const sci = formatScientific(rawVal);
                        pairs.push(`${k}: ${dec}   =   ${sci.str}`);
                    }else{
                        pairs.push(`${k}: ${String(rawVal)}`);
                    }
                });
                answerDiv.textContent = pairs.join('\n');
            }catch(e){ answerDiv.textContent = e.message || 'Error computing formula.'; console.error(e); }
        });

        const controls = document.createElement('div');
        controls.className = 'controls';
        controls.appendChild(computeBtn);
        paramsContainer.appendChild(controls);
        answerDiv.textContent = 'Enter values and press Compute.';
    });

    clearBtn.addEventListener('click', ()=>{
        paramsContainer.innerHTML = '';
        answerDiv.textContent = 'Result will appear here';
        formSel.selectedIndex = -1;
        catSel.selectedIndex = -1;
    });

    if(catSel.options.length>0){ catSel.selectedIndex = 0; populateFormulas(); }
    // --- Calculator integration ---
    const calcDisplay = document.getElementById('result'); // calculator display (id renamed)
    const calcButtons = document.querySelectorAll('#calc_buttons button');
    let focusedInput = null;
    // degree/radian mode for trig functions (false = radians)
    let degreeMode = false;
    const degToggleEl = document.getElementById('deg_toggle');
    if(degToggleEl){
        // default label; the button in HTML may show something else
        degToggleEl.textContent = 'RAD';
        degToggleEl.addEventListener('click', (e)=>{
            degreeMode = !degreeMode;
            degToggleEl.textContent = degreeMode ? 'DEG' : 'RAD';
            degToggleEl.classList.toggle('active', degreeMode);
        });
    }

    // set focusedInput when any parameter input receives focus
    paramsContainer.addEventListener('focusin', (e)=>{
        const t = e.target;
        if(t && t.tagName === 'INPUT' && t.dataset && t.dataset.key){
            focusedInput = t;
            // mirror current input into calc display for convenience
            try{ calcDisplay.textContent = t.value || ''; }catch(e){}
        }
    });

    // when inputs lose focus, keep the last focusedInput but don't clear it
    paramsContainer.addEventListener('focusout', (e)=>{
        // no-op; we keep focusedInput until another input is focused
    });

    function insertAtCursor(el, text){
        if(!el) return;
        const start = el.selectionStart || 0;
        const end = el.selectionEnd || 0;
        const value = el.value || '';
        const newVal = value.substring(0, start) + text + value.substring(end);
        el.value = newVal;
        const pos = start + text.length;
        el.setSelectionRange(pos, pos);
    }

    calcButtons.forEach(btn=>{
        btn.addEventListener('click', ()=>{
            const v = btn.dataset.value;
            if(!v) return;
            // handle DEG toggle button specially (don't insert text)
            if(v === 'deg' || btn.id === 'deg_toggle'){
                degreeMode = !degreeMode;
                try{ btn.textContent = degreeMode ? 'DEG' : 'RAD'; }catch(e){}
                btn.classList.toggle('active', degreeMode);
                if(degToggleEl && degToggleEl !== btn) degToggleEl.textContent = degreeMode ? 'DEG' : 'RAD';
                return;
            }
            if(v === 'AC'){
                // clear
                if(focusedInput){ focusedInput.value = ''; focusedInput.focus(); }
                calcDisplay.textContent = '';
                return;
            }

            if(v === '='){
                // evaluate expression either into focused input or display
                const exprSource = (focusedInput && focusedInput.value.trim() !== '') ? focusedInput.value : (calcDisplay.textContent || '');
                const result = (function(){
                    // try numeric parse first
                    const n1 = parseNumericInput(exprSource);
                    if(Number.isFinite(n1)) return n1;
                    const n2 = (function(){
                        try{ // reuse safe evaluator if available
                            if(typeof safeEvalExpression === 'function') return safeEvalExpression(exprSource);
                        }catch(e){}
                        return NaN;
                    })();
                    if(Number.isFinite(n2)) return n2;
                    return NaN;
                })();
                if(!Number.isFinite(result)){
                    calcDisplay.textContent = 'Err';
                }else{
                    // write back to focused input if present else update display
                    const pretty = (Math.abs(result) < 1e-6 || Math.abs(result) > 1e6) ? formatScientific(result).str : Number.parseFloat(result.toFixed(6)).toString();
                    if(focusedInput){
                        focusedInput.value = String(result);
                        calcDisplay.textContent = String(result);
                    }else{
                        calcDisplay.textContent = pretty;
                    }
                }
                return;
            }

            // for other buttons, append to focused input (if present) else to display
            if(focusedInput){
                // insert at cursor
                insertAtCursor(focusedInput, v);
                // mirror to display
                calcDisplay.textContent = focusedInput.value;
                focusedInput.focus();
            }else{
                calcDisplay.textContent = (calcDisplay.textContent || '') + v;
            }
        });
    });

    // keyboard shortcuts: Enter to evaluate, Escape to clear focused input / display
    document.addEventListener('keydown', (e)=>{
        if(e.key === 'Enter'){
            e.preventDefault();
            const exprSource = (focusedInput && focusedInput.value.trim() !== '') ? focusedInput.value : (calcDisplay.textContent || '');
            const n1 = parseNumericInput(exprSource);
            let result = NaN;
            if(Number.isFinite(n1)) result = n1;
            else result = safeEvalExpression(exprSource);
            if(!Number.isFinite(result)){
                calcDisplay.textContent = 'Err';
            }else{
                if(focusedInput){ focusedInput.value = String(result); calcDisplay.textContent = String(result); }
                else calcDisplay.textContent = (Math.abs(result) < 1e-6 || Math.abs(result) > 1e6) ? formatScientific(result).str : Number.parseFloat(result.toFixed(6)).toString();
            }
        }
        if(e.key === 'Escape'){
            if(focusedInput){ focusedInput.value = ''; focusedInput.focus(); }
            calcDisplay.textContent = '';
        }
    });
    // --- AI Assistant ---
    const aiToggle = document.getElementById('ai-toggle');
    const aiPanel = document.getElementById('ai-panel');
    const aiClose = document.getElementById('ai-close');
    const aiMessages = document.getElementById('ai-messages');
    const aiForm = document.getElementById('ai-form');
    const aiInput = document.getElementById('ai-input');
    const aiSettingsToggle = document.getElementById('ai-settings-toggle');
    const aiSettingsPanel = document.getElementById('ai-settings');
    const aiApiKeyInput = document.getElementById('ai-api-key');
    const aiSaveKey = document.getElementById('ai-save-key');
    const aiKeyStatus = document.getElementById('ai-key-status');

    // API key management
    let geminiApiKey = localStorage.getItem('gemini_api_key') || '';
    if (geminiApiKey) {
        aiApiKeyInput.value = geminiApiKey;
    }

    // Conversation history for Gemini context
    let conversationHistory = [];

    aiSettingsToggle.addEventListener('click', () => {
        aiSettingsPanel.classList.toggle('hidden');
        updateKeyStatus();
    });

    aiSaveKey.addEventListener('click', () => {
        const key = aiApiKeyInput.value.trim();
        if (key) {
            geminiApiKey = key;
            localStorage.setItem('gemini_api_key', key);
            aiKeyStatus.textContent = 'API key saved!';
            aiKeyStatus.className = 'ai-key-status success';
            conversationHistory = [];
        } else {
            geminiApiKey = '';
            localStorage.removeItem('gemini_api_key');
            aiKeyStatus.textContent = 'API key removed. Using offline mode.';
            aiKeyStatus.className = 'ai-key-status info';
            conversationHistory = [];
        }
    });

    function updateKeyStatus() {
        if (geminiApiKey) {
            aiKeyStatus.textContent = 'Key saved. Using Gemini AI.';
            aiKeyStatus.className = 'ai-key-status success';
        } else {
            aiKeyStatus.textContent = 'No key set. Using offline mode.';
            aiKeyStatus.className = 'ai-key-status info';
        }
    }

    aiToggle.addEventListener('click', () => {
        aiPanel.classList.toggle('hidden');
        if (!aiPanel.classList.contains('hidden')) {
            aiInput.focus();
            if (aiMessages.children.length === 0) {
                const mode = geminiApiKey ? 'Gemini AI' : 'offline mode';
                addAiMessage('assistant', `Hello! I can help you solve math and science problems (${mode}).\n\nTry asking something like:\n- "What is the area of a triangle with base 5 and height 10?"\n- "Explain Newton\'s second law"\n- "Convert 100 Fahrenheit to Celsius"\n- "sqrt(144) + 3^2"` + (!geminiApiKey ? '\n\nTip: Click the gear icon to add a free Gemini API key for smarter AI responses!' : ''));
            }
        }
    });
    aiClose.addEventListener('click', () => aiPanel.classList.add('hidden'));

    function addAiMessage(role, text) {
        const div = document.createElement('div');
        div.className = `ai-msg ${role}`;
        div.textContent = text;
        aiMessages.appendChild(div);
        aiMessages.scrollTop = aiMessages.scrollHeight;
        return div;
    }

    function addLoadingMessage() {
        const div = document.createElement('div');
        div.className = 'ai-msg loading';
        div.textContent = 'Thinking';
        div.id = 'ai-loading';
        aiMessages.appendChild(div);
        aiMessages.scrollTop = aiMessages.scrollHeight;
        return div;
    }

    function removeLoadingMessage() {
        const el = document.getElementById('ai-loading');
        if (el) el.remove();
    }

    // --- Gemini API integration ---
    const GEMINI_SYSTEM_PROMPT = `You are a helpful math and science assistant embedded in a Formula Solver web app. You help users solve mathematical problems, physics equations, chemistry questions, and unit conversions.

Available formula categories in the app: Triangle, Rectangle, Circle, Intersections, Electrical, Mechanics, Energy, Waves, Optics, Finance, Converter, Miscellaneous.

Guidelines:
- Give clear, step-by-step solutions
- Show the formula used and the substitution of values
- Provide the final numerical answer
- Keep responses concise but thorough
- Use plain text formatting (no markdown)
- For unit conversions, show the conversion factor used`;

    async function queryGemini(userMessage) {
        conversationHistory.push({ role: 'user', parts: [{ text: userMessage }] });

        const body = {
            system_instruction: { parts: [{ text: GEMINI_SYSTEM_PROMPT }] },
            contents: conversationHistory,
            generationConfig: {
                temperature: 0.7,
                maxOutputTokens: 1024,
            }
        };

        const url = `https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key=${geminiApiKey}`;

        const response = await fetch(url, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(body),
        });

        if (!response.ok) {
            const errData = await response.json().catch(() => ({}));
            const errMsg = errData?.error?.message || response.statusText;
            // Remove the failed user message from history
            conversationHistory.pop();
            throw new Error(errMsg);
        }

        const data = await response.json();
        const text = data?.candidates?.[0]?.content?.parts?.[0]?.text;
        if (!text) {
            conversationHistory.pop();
            throw new Error('Empty response from Gemini');
        }

        conversationHistory.push({ role: 'model', parts: [{ text }] });

        // Keep history manageable (last 20 messages)
        if (conversationHistory.length > 20) {
            conversationHistory = conversationHistory.slice(-20);
        }

        return text;
    }

    // --- Local offline fallback (pattern matching) ---
    const categoryKeywords = {
        triangle: ['triangle', 'pythagorean', 'hypotenuse', 'sine rule', 'cosine rule', 'base.*height'],
        rectangle: ['rectangle', 'rectangular', 'perimeter.*width', 'width.*height'],
        circle: ['circle', 'circumference', 'radius', 'sector', 'arc'],
        electrical: ['ohm', 'voltage', 'current', 'resistance', 'resistor', 'circuit', 'electrical', 'electric', 'power.*voltage', 'power.*current'],
        mechanics: ['force', 'momentum', 'work.*force', 'pressure', 'density', 'acceleration', 'velocity', 'displacement', 'hooke', 'spring', 'centripetal', 'newton'],
        energy: ['kinetic', 'potential energy', 'energy.*mass', 'energy.*velocity', 'energy.*height'],
        waves: ['wave', 'wavelength', 'frequency.*wavelength'],
        optics: ['lens', 'focal', 'optic', 'image distance', 'object distance'],
        finance: ['interest', 'principal', 'finance', 'loan', 'investment'],
        converter: ['convert', 'celsius', 'fahrenheit', 'kilometer', 'mile', 'degree.*radian', 'radian.*degree', 'kilogram', 'pound', 'liter', 'gallon', 'joule.*calorie', 'calorie.*joule', 'watt.*horsepower', 'horsepower.*watt', 'newton.*pound', 'cubic meter', 'cubic feet'],
        misc: ['frequency', 'period', 'density.*mass.*volume'],
        intersections: ['intersection', 'line.*line', 'line.*parabola', 'line.*circle', 'parabola.*parabola', 'circle.*circle'],
    };

    const formulaKeywords = {
        'triangle:area': ['area.*triangle', 'triangle.*area', 'base.*height'],
        'triangle:pythagorean': ['pythagorean', 'hypotenuse', 'pythagoras'],
        'triangle:sine rule': ['sine rule', 'law of sines', 'sin rule'],
        'triangle:cosine rule': ['cosine rule', 'law of cosines', 'cos rule'],
        'rectangle:area': ['area.*rectangle', 'rectangle.*area', 'width.*height.*area'],
        'rectangle:perimeter': ['perimeter.*rectangle', 'rectangle.*perimeter'],
        'circle:area': ['area.*circle', 'circle.*area', 'pi.*r.*squared', 'pi.*radius'],
        'circle:circumference': ['circumference', 'circle.*perimeter'],
        'circle:sector_area': ['sector', 'arc.*area'],
        'electrical:ohms_law': ['ohm', 'v.*i.*r', 'voltage.*current.*resistance'],
        'electrical:power': ['power.*voltage', 'power.*current', 'p.*v.*i'],
        'mechanics:force': ['force.*mass.*acceleration', 'f.*m.*a', 'newton.*second'],
        'mechanics:momentum': ['momentum'],
        'mechanics:work': ['work.*force.*distance'],
        'mechanics:pressure': ['pressure.*force.*area'],
        'mechanics:density': ['density.*mass.*volume'],
        'mechanics:acceleration': ['acceleration.*velocity.*time', 'delta.*v.*delta.*t'],
        'mechanics:velocity': ['velocity.*distance.*time', 'speed.*distance.*time'],
        'mechanics:displacement': ['displacement.*initial.*velocity', 's.*u.*a.*t'],
        'mechanics:hookes_law': ['hooke', 'spring.*constant', 'spring.*force'],
        'mechanics:centripetal_force': ['centripetal'],
        'energy:kinetic': ['kinetic energy', 'kinetic', 'ke.*mass.*velocity'],
        'energy:potential': ['potential energy', 'pe.*mass.*height', 'mgh'],
        'waves:wave_speed': ['wave speed', 'wave.*frequency.*wavelength'],
        'optics:lens': ['lens', 'focal length', 'image distance', 'object distance'],
        'finance:simple_interest': ['simple interest', 'interest.*principal.*rate.*time'],
        'converter:celsius_fahrenheit': ['celsius.*fahrenheit', 'fahrenheit.*celsius', 'temperature.*convert'],
        'converter:km_miles': ['kilometer.*mile', 'mile.*kilometer', 'km.*mile'],
        'converter:degrees_radians': ['degree.*radian', 'radian.*degree'],
        'converter:kg_pounds': ['kilogram.*pound', 'pound.*kilogram', 'kg.*lb', 'lb.*kg'],
        'converter:liters_gallons': ['liter.*gallon', 'gallon.*liter'],
        'converter:newtons_lbf': ['newton.*pound', 'pound.*force.*newton'],
        'converter:watts_hp': ['watt.*horsepower', 'horsepower.*watt', 'hp.*watt'],
        'converter:joules_calories': ['joule.*calorie', 'calorie.*joule'],
    };

    function extractNumbers(text) {
        const numbers = [];
        const pattern = /(?:^|[\s=:,])([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)/g;
        let match;
        while ((match = pattern.exec(text)) !== null) {
            numbers.push(parseFloat(match[1]));
        }
        return numbers;
    }

    function extractNamedValues(text) {
        const result = {};
        const patterns = [
            /(\w[\w\s]*?)\s*(?:is|=|equals|of)\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)/gi,
            /(\w[\w\s]*?)\s*:\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)/gi,
        ];
        for (const pat of patterns) {
            let m;
            while ((m = pat.exec(text)) !== null) {
                result[m[1].trim().toLowerCase()] = parseFloat(m[2]);
            }
        }
        return result;
    }

    function matchFormula(text) {
        const lower = text.toLowerCase();
        let bestMatch = null;
        let bestScore = 0;

        for (const [key, keywords] of Object.entries(formulaKeywords)) {
            for (const kw of keywords) {
                const regex = new RegExp(kw, 'i');
                if (regex.test(lower)) {
                    const [cat, fid] = key.split(':');
                    const formula = (formulaCatalog[cat] || []).find(f => f.id === fid);
                    if (formula) {
                        const score = kw.length + 10;
                        if (score > bestScore) {
                            bestScore = score;
                            bestMatch = { category: cat, formula };
                        }
                    }
                }
            }
        }

        if (bestMatch) return bestMatch;

        for (const [cat, keywords] of Object.entries(categoryKeywords)) {
            for (const kw of keywords) {
                const regex = new RegExp(kw, 'i');
                if (regex.test(lower)) {
                    const formulas = formulaCatalog[cat];
                    if (formulas && formulas.length > 0) {
                        return { category: cat, formula: formulas[0] };
                    }
                }
            }
        }

        return null;
    }

    function mapValuesToParams(formula, text) {
        const named = extractNamedValues(text);
        const vals = {};

        for (const param of formula.params) {
            const key = param.key;
            const label = param.label.toLowerCase();

            for (const [name, value] of Object.entries(named)) {
                if (label.includes(name) || key.toLowerCase() === name || name.includes(key.toLowerCase())) {
                    vals[key] = value;
                    break;
                }
            }
        }

        if (Object.keys(vals).length === 0) {
            const numbers = extractNumbers(text);
            const fillableParams = formula.params.filter(p => p.type === 'number');
            for (let i = 0; i < Math.min(numbers.length, fillableParams.length); i++) {
                vals[fillableParams[i].key] = numbers[i];
            }
        }

        return vals;
    }

    function processOfflineQuery(text) {
        const lower = text.toLowerCase().trim();

        if (/^(hi|hello|hey|greetings|good\s*(morning|afternoon|evening))[\s!.]*$/i.test(lower)) {
            return 'Hello! I\'m in offline mode. I can solve formulas and evaluate expressions. Add a Gemini API key (gear icon) for smarter AI responses!';
        }

        if (/^(help|what can you do|how do you work)[\s?]*$/i.test(lower)) {
            return 'Offline mode can:\n1. Solve formulas (triangles, physics, etc.)\n2. Evaluate math expressions\n3. Calculate percentages\n\nFor full AI capabilities, add a free Gemini API key via the gear icon!';
        }

        // Percentage
        let m = lower.match(/(\d+\.?\d*)\s*(%\s*(?:of)|percent\s*(?:of))\s*(\d+\.?\d*)/);
        if (m) {
            const pct = parseFloat(m[1]);
            const total = parseFloat(m[3]);
            return `${pct}% of ${total} = ${(pct / 100) * total}`;
        }

        // Formula matching
        const match = matchFormula(text);
        if (match) {
            const { category, formula } = match;
            const vals = mapValuesToParams(formula, text);
            const filledCount = Object.values(vals).filter(v => Number.isFinite(v)).length;

            if (filledCount === 0) {
                return `Formula: ${formula.name}\nCategory: ${category}\n\nProvide values like: ${formula.params.map(p => p.label.split('(')[0].trim()).join(', ')}`;
            }

            try {
                const solved = formula.solve(vals);
                const lines = [`Formula: ${formula.name}`, `Category: ${category}`, '', 'Given:'];
                for (const param of formula.params) {
                    if (vals[param.key] !== undefined && Number.isFinite(vals[param.key])) {
                        lines.push(`  ${param.label} = ${vals[param.key]}`);
                    }
                }
                lines.push('', 'Result:');
                for (const [key, value] of Object.entries(solved)) {
                    if (typeof value === 'number') {
                        const paramDef = formula.params.find(p => p.key === key);
                        const label = paramDef ? paramDef.label : key;
                        const dec = Number.parseFloat(value.toFixed(6));
                        lines.push(`  ${label} = ${dec}`);
                    }
                }
                return lines.join('\n');
            } catch (e) {
                return `Formula: ${formula.name}\nError: ${e.message}`;
            }
        }

        // Direct expression
        const cleaned = text.replace(/^(?:what is|calculate|compute|evaluate|solve|find)\s*/i, '').trim();
        if (cleaned) {
            const result = safeEvalExpression(cleaned);
            if (Number.isFinite(result)) {
                return `${cleaned} = ${result}`;
            }
        }

        return 'I couldn\'t understand that (offline mode). Try:\n- "Area of triangle with base 5 and height 10"\n- "sqrt(144) + 3^2"\n- "15% of 200"\n\nFor better results, add a free Gemini API key!';
    }

    // --- Submit handler ---
    let isProcessing = false;

    aiForm.addEventListener('submit', async (e) => {
        e.preventDefault();
        const text = aiInput.value.trim();
        if (!text || isProcessing) return;

        addAiMessage('user', text);
        aiInput.value = '';

        if (geminiApiKey) {
            // Use Gemini API
            isProcessing = true;
            aiInput.disabled = true;
            const loadingEl = addLoadingMessage();

            try {
                const response = await queryGemini(text);
                removeLoadingMessage();
                addAiMessage('assistant', response);
            } catch (err) {
                removeLoadingMessage();
                const errMsg = err.message || 'Unknown error';
                if (errMsg.includes('API key') || errMsg.includes('api key') || errMsg.includes('403') || errMsg.includes('401')) {
                    addAiMessage('assistant', `API key error: ${errMsg}\n\nFalling back to offline mode for this question.`);
                    addAiMessage('assistant', processOfflineQuery(text));
                } else {
                    addAiMessage('assistant', `API error: ${errMsg}\n\nUsing offline fallback:`);
                    addAiMessage('assistant', processOfflineQuery(text));
                }
            } finally {
                isProcessing = false;
                aiInput.disabled = false;
                aiInput.focus();
            }
        } else {
            // Offline mode
            setTimeout(() => {
                addAiMessage('assistant', processOfflineQuery(text));
            }, 100);
        }
    });

    // Prevent Enter/Escape in AI input from triggering global handlers
    aiInput.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') e.stopPropagation();
        if (e.key === 'Escape') {
            e.stopPropagation();
            aiPanel.classList.add('hidden');
        }
    });
    // Also stop propagation from the API key input
    aiApiKeyInput.addEventListener('keydown', (e) => {
        e.stopPropagation();
    });
});
console.log('Formula Solver script loaded.');
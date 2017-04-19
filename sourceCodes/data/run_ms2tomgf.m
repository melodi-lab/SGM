for id = 1:100
ms2file = ['ms2data/human/Linfeng-', int2str(id) ,'.ms2'];
for charge = 1:5
    mgfFile = ['mascotdata/human/Linfeng-', int2str(id) ,'-ch',int2str(charge),'.mgf'];
    charge_filter= charge;
    numPeps = 2000;
    ms2tomgf(ms2file, mgfFile, charge_filter, numPeps);
end 
end

for id = 1:20
ms2file = ['ms2data/plasm/plasm-', int2str(id) ,'.ms2'];
for charge = 1:5
    mgfFile = ['mascotdata/plasm/plasm-', int2str(id) ,'-ch',int2str(charge),'.mgf'];
    charge_filter= charge;
    numPeps = 2000;
    ms2tomgf(ms2file, mgfFile, charge_filter, numPeps);
end 
end

for id = 1:3
ms2file = ['ms2data/yeast/yeast-0', int2str(id) ,'.ms2'];
for charge = 1:5
    mgfFile = ['mascotdata/yeast/yeast-', int2str(id) ,'-ch',int2str(charge),'.mgf'];
    charge_filter= charge;
    numPeps = 2000;
    ms2tomgf(ms2file, mgfFile, charge_filter, numPeps);
end 
end

for id = 1:3
ms2file = ['ms2data/worm/worm-0', int2str(id) ,'.ms2'];
for charge = 1:5
    mgfFile = ['mascotdata/worm/worm-', int2str(id) ,'-ch',int2str(charge),'.mgf'];
    charge_filter= charge;
    numPeps = 2000;
    ms2tomgf(ms2file, mgfFile, charge_filter, numPeps);
end 
end





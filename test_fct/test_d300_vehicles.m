
DrugNames = {
    'DMSO_based'     ''
    'Brij35_based'     ''
    'Brij35+gly_based'     ''
    'TrytonX100_based'     ''
    'TrytonX100+gly_based'     ''
    'Tween20_based'     ''
    'Tween20+gly_based'     ''
    };
Vehicles = {'DMSO' 'brij35' 'brij35+glycerol' 'tryton' 'tryton+glycerol' 'tween' 'tween+glycerol'};

StockConc = 5e3*ones(1,size(DrugNames,1));
SingleDoses = repmat({10.^[-3:.5:1]},size(DrugNames,1),1);

Design = TreatmentDesign(DrugNames(:,1), DrugNames(:,2), ...
    SingleDoses, 3, 'Seed', 1, 'stock_conc', StockConc, ...
    'Vehicle', Vehicles);

t_plates = table({'plate1';'plate2'}, repmat({'test.mat'},2,1), [1;2], ...
    'variablenames', {'Barcode' 'TreatmentFile' 'DesignNumber'});

hpdd_exporter('test_D300', Design, t_plates);
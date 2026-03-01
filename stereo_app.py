def render_smart_2d(mol):
    # نسخة من الموليكيول عشان منبوظش الأصل
    m = Chem.Mol(mol)
    is_allene = m.HasSubstructMatch(Chem.MolFromSmarts("C=C=C"))
    
    # 1. ضروري جداً إضافة الهيدروجين للألينات عشان تترسم الـ Wedges صح
    m = Chem.AddHs(m)
    
    # 2. توليد إحداثيات 3D الأول عشان نعرف الـ Stereochemistry
    if AllChem.EmbedMolecule(m, AllChem.ETKDG()) != -1:
        # 3. تحويل إحداثيات الـ 3D لـ 2D مع الحفاظ على الـ Stereochemistry
        AllChem.Compute2DCoords(m)
        # 4. تحديد روابط الـ Wedge بناءً على الـ Conformer
        Chem.WedgeMolBonds(m, m.GetConformer())
    else:
        # لو فشل الـ Embedding ارسم عادي
        AllChem.Compute2DCoords(m)

    d_opts = Draw.MolDrawOptions()
    d_opts.addStereoAnnotation = True
    
    # 5. تحسينات شكلية للألين
    if is_allene:
        d_opts.bondLineWidth = 3.0    # خطوط سميكة عشان الـ Wedges تبان
        d_opts.minFontSize = 18
    else:
        d_opts.bondLineWidth = 1.6
        d_opts.minFontSize = 14

    img = Draw.MolToImage(m, size=(500, 500), options=d_opts)
    return img

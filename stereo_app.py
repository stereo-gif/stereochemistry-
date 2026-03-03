from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def render_smart_2d(mol):
    # 1. نسخة من الموليكيول عشان الأصل ميتأثرش
    m = Chem.Mol(mol)
    
    # التأكد من وجود ألين (C=C=C) لضبط التنسيق لاحقاً
    is_allene = m.HasSubstructMatch(Chem.MolFromSmarts("C=C=C"))
    
    # 2. إضافة الهيدروجين ضروري جداً للألينات لظهور الـ Wedges بشكل صحيح
    m = Chem.AddHs(m)
    
    # 3. محاولة توليد إحداثيات 3D لمعرفة الـ Stereochemistry (خاصة الـ Axial)
    if AllChem.EmbedMolecule(m, AllChem.ETKDG()) != -1:
        # تحويل الـ 3D لـ 2D مع الحفاظ على الاتجاهات الفراغية
        AllChem.Compute2DCoords(m)
        # رسم الروابط بنظام المثلثات (Wedge/Dash) بناءً على الـ Conformer
        Chem.WedgeMolBonds(m, m.GetConformer())
    else:
        # لو فشل الـ 3D (زي في الجزيئات العملاقة) نرسم 2D عادي
        AllChem.Compute2DCoords(m)

    # 4. إعدادات الرسم (Options)
    d_opts = Draw.MolDrawOptions()
    
    # --- حل مشكلة Ra / Sa و isomer 2 ---
    d_opts.addStereoAnnotation = False  # ده اللي بيلغي الكتابة النصية فوق الذرات
    
    # ضبط الخطوط بناءً على نوع المركب
    if is_allene:
        d_opts.bondLineWidth = 3.0    # خطوط سميكة للوضوح
        d_opts.minFontSize = 18
    else:
        d_opts.bondLineWidth = 1.6
        d_opts.minFontSize = 14

    # 5. توليد الصورة النهائية
    # الـ legend="" بتضمن إن مفيش أي اسم (زي Search Bitter) يظهر تحت أو فوق المركب
    img = Draw.MolToImage(m, size=(500, 500), options=d_opts, legend="")
    
    return img

# --- مثال لتشغيل الكود ---
# سنستخدم ألين (Allene) كمثال لأنه الأصعب في الرسم
smiles_string = "CC=C=CC" 
test_mol = Chem.MolFromSmiles(smiles_string)

if test_mol:
    result_img = render_smart_2d(test_mol)
    result_img.show() # هيفتح الصورة عندك على الجهاز

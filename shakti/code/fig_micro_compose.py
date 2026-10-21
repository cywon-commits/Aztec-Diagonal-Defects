"""Fig. 16: stack the three process panels drawn by draw.py (process_figures/*.png)."""
import os
from PIL import Image, ImageDraw, ImageFont
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(HERE, '..', 'process_figures')
ims = [Image.open(os.path.join(SRC, f)).convert('RGB') for f in ['single_pair_creation.png', 'defect_hop.png', 'two_pair_creation.png']]
ims = [im.crop((0, 0, im.width, im.height - 55)) for im in ims[:2]] + [ims[2]]   # keep the legend line only once
W = max(i.width for i in ims); out = Image.new('RGB', (W, sum(i.height for i in ims) + 30), 'white'); d = ImageDraw.Draw(out)
try:
    font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 34)
except OSError:
    font = None
y = 0
for k, im in enumerate(ims):
    out.paste(im, (0, y)); d.text((12, y + 8), '(%s)' % 'abc'[k], fill='black', font=font); y += im.height + 10
out.save(os.path.join(ROOT, 'figures', 'fig16_micro.png')); print('wrote figures/fig16_micro.png')

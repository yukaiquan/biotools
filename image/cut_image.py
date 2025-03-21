from PIL import Image
import sys
import os

def combine_images(image_paths, rows, cols, output_path):
    # 计算每张图片的宽度和高度
    widths, heights = [], []
    for path in image_paths:
        with Image.open(path) as img:
            widths.append(img.width)
            heights.append(img.height)
    
    # 确定最大宽度和高度
    max_width = max(widths)
    max_height = max(heights)
    
    # 计算合并图片的总宽度和高度
    total_width = cols * max_width
    total_height = rows * max_height
    
    # 创建新画布
    combined_image = Image.new('RGB', (total_width, total_height))
    
    # 粘贴图片
    for idx, path in enumerate(image_paths):
        with Image.open(path) as img:
            # 计算每张图片的放置位置
            col = idx % cols
            row = idx // cols
            # 调整图片大小并粘贴到画布上
            combined_image.paste(img.resize((max_width, max_height)), (col * max_width, row * max_height))
    
    # 保存合并后的图片
    combined_image.save(output_path)

def crop_image(image_path, crop_area, output_path,angle):
    with Image.open(image_path) as img:
        rotated_image = img.rotate(-angle, expand=True)
        # 裁剪图片
        cropped_image = rotated_image.crop(crop_area)
        # 逆时针90度旋转
        cropped_image = cropped_image.transpose(Image.ROTATE_90)
        # 保存裁剪后的图片
        cropped_image.save(output_path)

# 图片路径
# image_path = 'path/to/your/image.jpg'
samples = sys.argv[1]
samples_list = []
with open(samples,'r') as f:
    for line in f.readlines():
        if line.strip() != '':
            samples_list.append(line.strip())

output_dir = sys.argv[2]
# 裁剪区域：(left, upper, right, lower)
crop_area = (913, 324, 2456, 2532)  # 举例区域
angle = 0.5

output_list = []
for sample in samples_list:
    name = sample.split('.')[0]
    output_name = os.path.join(output_dir, name + '_cropped.bmp')
    if os.path.exists(output_name):
        print('skip',output_name)
    else:
        crop_image(sample, crop_area, output_name,angle)
    output_list.append(output_name)

# # 行数和列数
# rows = 10
# cols = 10
# # 输出路径
# output_combine_path = sys.argv[3]
# combine_images(output_list, rows, cols, output_combine_path)


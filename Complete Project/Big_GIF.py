from PIL import Image, ImageDraw
import imageio.v2 as imageio
import os
'''
# Generate list of filenames
steps = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
gif_paths = [f"pressure_vs_density_rinf{i}.gif" for i in steps]

# Load each GIF using get_reader
gif_readers = [imageio.get_reader(gif) for gif in gif_paths]

# Get size from the first frame of the first GIF
first_frame = gif_readers[0].get_data(0)
height, width, _ = first_frame.shape  # shape is (height, width, channels)

# Grid layout: 5 columns x 2 rows
cols, rows = 5, 2
final_width = cols * width
final_height = rows * height

# Line thickness and color
line_color = (0, 0, 0)  # Black line
line_thickness = 2  # pixels

# Prepare output frames
output_frames = []

# Find the max length among all GIFs
max_length = max(reader.get_length() for reader in gif_readers)

for frame_idx in range(max_length):
    grid_img = Image.new('RGB', (final_width, final_height))
    draw = ImageDraw.Draw(grid_img)

    for idx, reader in enumerate(gif_readers):
        try:
            frame = reader.get_data(frame_idx)
        except (IndexError, EOFError):
            # If out of range, repeat the last frame
            frame = reader.get_data(reader.get_length() - 1)

        img = Image.fromarray(frame)
        x = (idx % cols) * width
        y = (idx // cols) * height
        grid_img.paste(img, (x, y))

    # Draw horizontal lines between rows
    for row in range(1, rows):
        y_line = row * height
        draw.line(
            [(0, y_line), (final_width, y_line)],
            fill=line_color,
            width=line_thickness
        )

    # Optional: Draw vertical lines between columns
    for col in range(1, cols):
        x_line = col * width
        draw.line([(x_line, 0), (x_line, final_height)], fill=line_color, width=line_thickness)

    output_frames.append(grid_img.convert('P', palette=Image.ADAPTIVE))

# Save the final combined GIF
output_frames[0].save(
    'combined_pressure_vs_density.gif',
    save_all=True,
    append_images=output_frames[1:],
    duration=300,
    loop=0,
    disposal=2
)
'''


# List of PNG filenames you want to combine
steps = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
png_paths = [f"ADM_mass_vs_omega_rinf_{i}.png" for i in steps]

# Load all images
images = [Image.open(png) for png in png_paths]

# Get image size from the first image (assuming all have same size)
width, height = images[0].size

# Grid layout: 5 columns x 2 rows
cols, rows = 5, 2
final_width = cols * width
final_height = rows * height

# Create a new blank image to paste all images into
combined_img = Image.new('RGB', (final_width, final_height))

# Paste images into the grid
for idx, img in enumerate(images):
    x = (idx % cols) * width
    y = (idx // cols) * height
    combined_img.paste(img, (x, y))

# Draw lines between images
draw = ImageDraw.Draw(combined_img)
line_color = (0, 0, 0)  # black
line_thickness = 2

# Horizontal lines
for row in range(1, rows):
    y_line = row * height
    draw.line([(0, y_line), (final_width, y_line)], fill=line_color, width=line_thickness)

# Vertical lines
for col in range(1, cols):
    x_line = col * width
    draw.line([(x_line, 0), (x_line, final_height)], fill=line_color, width=line_thickness)

# Save combined image
combined_img.save('combined_ADM_mass_vs_omega_grid.png')


print("✅")
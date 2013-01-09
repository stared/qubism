// 100914, easyim
#include"easyim.h"

void EI_Start()
{
     imlib_set_cache_size(2048 * 1024);
     imlib_context_set_display(EXInfo.display);
     imlib_context_set_visual(EXInfo.visual);
     imlib_context_set_colormap(EXInfo.colormap);
     if (EXCW->use_buffer)
	  imlib_context_set_drawable(EXCW->buffer); 
     else
	  imlib_context_set_drawable(EXCW->window); 
}

EI_Image EI_Load(const char *name) // also sets current image!
{
     EI_Image image;
     image = imlib_load_image(name);
     imlib_context_set_image(image);
     return image;
}

int EI_Get_Width()
{
     return imlib_image_get_width();
}

int EI_Get_Height()
{
     return imlib_image_get_height();
}

void EI_Render(int x, int y)
{
     imlib_render_image_on_drawable(x,y);
}

void EI_Render_Scaled(int x, int y, int w, int h)
{
     EI_Image nueva=imlib_create_image(w,h);
     EI_Image antigua=imlib_context_get_image();
     int oldw=EI_Get_Width();
     int oldh=EI_Get_Height();
     imlib_context_set_image(nueva);
     imlib_blend_image_onto_image(antigua,0,0,0,oldw,oldh,0,0,w,h);
     imlib_render_image_on_drawable(x,y);
     imlib_free_image();
     imlib_context_set_image(antigua);
}

void EI_Free()
{
    imlib_free_image();
}

EI_Image EI_Capture(int x, int y, int w, int h) // also sets current image
{
     EI_Image image=imlib_create_image_from_drawable(0,0,0,w,h,1);
     imlib_context_set_image(image);
     return image;
}

void EI_Save(const char *name)
{
     imlib_save_image(name);
}



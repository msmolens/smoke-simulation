#include "ImageWriter.h"

ImageWriter::ImageWriter(int _w, int _h) : w(_w), h(_h) { }

void
ImageWriter::set_dimensions(int _w, int _h)
{
	w = _w;
	h = _h;
}

void
ImageWriter::output (int frameNum)
{
//	int w=640;
//	int h=480;

	unsigned char
		array[h*w*3];

	glReadPixels ( 0,
		       0,
		       w,
		       h,
		       GL_RGB,
		       GL_UNSIGNED_BYTE,
		       (GLvoid*)array );

	ImageInfo* imageInfo;
	Image* image;
	ExceptionInfo exceptionInfo;


	GetExceptionInfo ( &exceptionInfo );
	imageInfo = CloneImageInfo ((ImageInfo*) NULL );
	snprintf ( imageInfo->filename, 20, "movie%04d.jpg", frameNum );
 
	image = ConstituteImage ( w, 
				  h,
				  "RGB",
				  CharPixel,
				  array,
				  &exceptionInfo );

	Image *f = FlipImage(image, &exceptionInfo);
	DestroyImage(image);

	snprintf ( f->filename, 20, "movie%04d.jpg", frameNum );

  
//	sprintf ( f->magick, "PNG" );
//	sprintf ( imageInfo->magick, "PNG" );
  
	if ( !WriteImage ( imageInfo, f ) ) {
		cerr << "Error writing image: " << image->exception.reason << endl;
	}
  
        //DestroyConstitute();
	DestroyImage ( f );
	DestroyImageInfo ( imageInfo );
	DestroyExceptionInfo(&exceptionInfo);
//        DestroyMagick();

}

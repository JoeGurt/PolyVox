#include "OpenGLWidget.h"

#include <QMouseEvent>
#include <QMatrix4x4>
#include <QCoreApplication>
#include <QTimer>
//#include <QtMath>

using namespace PolyVox;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Public functions
////////////////////////////////////////////////////////////////////////////////
OpenGLWidget::OpenGLWidget(QWidget *parent)
	:QGLWidget(parent)
{
}

void OpenGLWidget::setCameraTransform(QVector3D position, float pitch, float yaw)
{
	mCameraPosition = position;
	mCameraYaw = yaw;
	mCameraPitch = pitch;
}

void OpenGLWidget::mousePressEvent(QMouseEvent* event)
{
	// Initialise these variables which will be used when the mouse actually moves.
	m_CurrentMousePos = event->pos();
	m_LastFrameMousePos = m_CurrentMousePos;
}

void OpenGLWidget::mouseMoveEvent(QMouseEvent* event)
{
	// Update the x and y rotations based on the mouse movement.
	m_CurrentMousePos = event->pos();
	QPoint diff = m_CurrentMousePos - m_LastFrameMousePos;
	mCameraYaw -= diff.x() * mCameraRotateSpeed;
	mCameraPitch -= diff.y() * mCameraRotateSpeed;
	m_LastFrameMousePos = m_CurrentMousePos;
}

void OpenGLWidget::keyPressEvent(QKeyEvent* event)
{
	if (event->key() == Qt::Key_Escape)
	{
		close();
	}

	mPressedKeys.append(event->key());
}

void OpenGLWidget::keyReleaseEvent(QKeyEvent* event)
{
	mPressedKeys.removeAll(event->key());
}

////////////////////////////////////////////////////////////////////////////////
// Protected functions
////////////////////////////////////////////////////////////////////////////////
void OpenGLWidget::initializeGL()
{
	if (!initializeOpenGLFunctions())
	{
		std::cerr << "Could not initialize OpenGL functions" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//Print out some information about the OpenGL implementation.
	std::cout << "OpenGL Implementation Details:" << std::endl;
	if(glGetString(GL_VENDOR))
	  std::cout << "\tGL_VENDOR: " << glGetString(GL_VENDOR) << std::endl;
	if(glGetString(GL_RENDERER))
	  std::cout << "\tGL_RENDERER: " << glGetString(GL_RENDERER) << std::endl;
	if(glGetString(GL_VERSION))
	  std::cout << "\tGL_VERSION: " << glGetString(GL_VERSION) << std::endl;
	if(glGetString(GL_SHADING_LANGUAGE_VERSION))
	  std::cout << "\tGL_SHADING_LANGUAGE_VERSION: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

	//Set up the clear colour
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearDepth(1.0f);
	
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glDepthFunc(GL_LEQUAL);
	glDepthRange(0.0, 1.0);

	initialize();

	// Start a timer to drive the main rendering loop.
	QTimer* timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this, SLOT(update()));
	timer->start(0);

	mElapsedTimer.start();
}

void OpenGLWidget::resizeGL(int w, int h)
{
	//Setup the viewport
	glViewport(0, 0, w, h);
	
	auto aspectRatio = w / (float)h;
	float zNear = 1.0;
	float zFar = 1000.0;
	
	projectionMatrix.setToIdentity();
	//projectionMatrix.frustum(-aspectRatio, aspectRatio, -1, 1, zNear, zFar);
	projectionMatrix.perspective(mCameraFOV, aspectRatio, zNear, zFar);
}

void OpenGLWidget::paintGL()
{
	// Direction : Spherical coordinates to Cartesian coordinates conversion
	QVector3D cameraForward(
		cos(mCameraPitch) * sin(mCameraYaw),
		sin(mCameraPitch),
		cos(mCameraPitch) * cos(mCameraYaw)
		);

	// Right vector
	QVector3D cameraRight(
		sin(mCameraYaw - 3.14f / 2.0f),
		0,
		cos(mCameraYaw - 3.14f / 2.0f)
		);

	// Up vector
	QVector3D cameraUp = QVector3D::crossProduct(cameraRight, cameraForward);

	// Get the elapsed time since last frame and convert to seconds.
	float deltaTime = mElapsedTimer.restart() / 1000.0f;

	// Move forward
	if ((mPressedKeys.contains(Qt::Key_Up)) || (mPressedKeys.contains(Qt::Key_W)))
	{
		mCameraPosition += cameraForward * deltaTime * mCameraMoveSpeed;
	}
	// Move backward
	if ((mPressedKeys.contains(Qt::Key_Down)) || (mPressedKeys.contains(Qt::Key_S)))
	{
		mCameraPosition -= cameraForward * deltaTime * mCameraMoveSpeed;
	}
	// Strafe right
	if ((mPressedKeys.contains(Qt::Key_Right)) || (mPressedKeys.contains(Qt::Key_D)))
	{
		mCameraPosition += cameraRight * deltaTime * mCameraMoveSpeed;
	}
	// Strafe left
	if ((mPressedKeys.contains(Qt::Key_Left)) || (mPressedKeys.contains(Qt::Key_A)))
	{
		mCameraPosition -= cameraRight * deltaTime * mCameraMoveSpeed;
	}

	viewMatrix.setToIdentity();
	viewMatrix.lookAt(
		mCameraPosition,           // Camera is here
		mCameraPosition + cameraForward, // and looks here : at the same position, plus "direction"
		cameraUp                  // Head is up (set to 0,-1,0 to look upside-down)
		);

	//Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	renderOneFrame();
	
	// Check for errors.
	GLenum errCode = glGetError();
	if(errCode != GL_NO_ERROR)
	{
	  std::cerr << "OpenGL Error: " << errCode << std::endl;
	}
}

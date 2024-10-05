import vtk
import sys
import os

def view(filename):
    colors = vtk.vtkNamedColors()
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()

    extractEdges = vtk.vtkExtractEdges()
    extractEdges.SetInputConnection(reader.GetOutputPort())

    edgeMapper = vtk.vtkPolyDataMapper()
    edgeMapper.SetInputConnection(extractEdges.GetOutputPort())
    edgeMapper.SetScalarRange(0, 1)

    edgeActor = vtk.vtkActor()
    edgeActor.SetMapper(edgeMapper)
    edgeActor.GetProperty().SetLineWidth(2)

    labelMapper = vtk.vtkLabeledDataMapper()
    labelMapper.SetInputConnection(reader.GetOutputPort())
    labelActor = vtk.vtkActor2D()
    labelActor.SetMapper(labelMapper)
    labelActor.GetProperty().SetPointSize(25)

    # extractPoints = vtk.vtkExtractPoints()
    # extractPoints.SetInputConnection(reader.GetOutputPort())

    contextView = vtk.vtkContextView()
    renderer = contextView.GetRenderer()
    renderWindow = contextView.GetRenderWindow()
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    renderer.AddActor(edgeActor)
    renderer.AddActor(labelActor)
    renderer.SetBackground(colors.GetColor3d('SlateGray'))

    axes = vtk.vtkAxesActor()

    widget = vtk.vtkOrientationMarkerWidget()
    rgba = [1] * 4
    colors.GetColor('Carrot', rgba)
    widget.SetOutlineColor(rgba[0], rgba[1], rgba[2])
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(renderWindowInteractor)
    widget.SetViewport(0.0, 0.0, 0.4, 0.4)
    widget.SetEnabled(1)
    widget.InteractiveOn()

    aCamera = vtk.vtkCamera()
    aCamera.Azimuth(-40.0)
    aCamera.Elevation(50.0)

    renderer.SetActiveCamera(aCamera)
    renderer.ResetCamera()

    renderWindow.SetSize(640, 480)
    renderWindow.SetWindowName('ReadLegacyUnstructuredGrid')
    renderWindow.Render()

    renderWindowInteractor.Start()

def main() :
    if len(sys.argv) < 2 :
        print('Error: please specify vtk file')
        return

    filename = sys.argv[1]    
    view(filename)

if __name__ == '__main__' :
    main()
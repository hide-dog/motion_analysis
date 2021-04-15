using Printf

function main()
    output_dir_name  = "post_result"
    input_result_dir = "result"

    #vtk_頭の文字
    front = "# vtk DataFile Version 2.0"*"\n"
    back  = "\n"*"ASCII"*"\n"*"DATASET UNSTRUCTURED_GRID"*"\n"
    write_file(output_dir_name, input_result_dir, front, back)
end
    
function write_file(outdir, inresult_dir, front, back)
    x = zeros(8,3)

    inf = readdir(inresult_dir)
    for i in 1:length(inf)
        if occursin(".dat", inf[i]) == true
            dname = replace(inf[i],".dat" => "")
            out_file = replace(inf[i],".dat" => ".vtk")
            print("start writing "*out_file*"\n")

            fff=[]
            open(inresult_dir*"/"*inf[i], "r") do f
                fff=read(f,String)
            end 
            fff = split(fff,"\n",keepempty=false)

            for j in 2:length(fff)
                fff[j] = replace(fff[j]," \r" => "")
                temp = split(fff[j]," ")
                x[j-1,1] = parse(Float64,temp[1])
                x[j-1,2] = parse(Float64,temp[2])
                x[j-1,3] = parse(Float64,temp[3])
            end
            
            write_points(x, out_file, outdir, dname, front, back)
            write_cells(out_file, outdir)
            write_result(out_file, outdir)
            print("fin writing "*out_file*"\n")
        end
    end
end 

function  write_points(x_arr, out_file, outdir, dname, front, back)
    a = size(x_arr)[1]
    a_st = @sprintf("%1.0f", a)

    fff = outdir*"/"*out_file
    open(fff,"w") do f
        write(f,front)
        write(f,dname)
        write(f,back)
        write(f,"POINTS "*a_st*" float\n")
        for i in 1:a
            x = @sprintf("%7.7f", x_arr[i,1])
            y = @sprintf("%7.7f", x_arr[i,2])
            z = @sprintf("%7.7f", x_arr[i,3])

            write(f, x*" "*y*" "*z*"\n")
        end
    end
    println("fin writing points")
end

function write_cells(out_file,outdir)
    
    fff = outdir*"/"*out_file
    open(fff,"a") do f
        write(f, "CELLS "*"1"*" "*"9"*"\n")
        
        d1 = @sprintf("%1.0f", 0)
        d2 = @sprintf("%1.0f", 1)
        d3 = @sprintf("%1.0f", 3)
        d4 = @sprintf("%1.0f", 2)
        d5 = @sprintf("%1.0f", 4)
        d6 = @sprintf("%1.0f", 5)
        d7 = @sprintf("%1.0f", 7)
        d8 = @sprintf("%1.0f", 6)
        write(f, "8 "*d1*" "*d2*" "*d3*" "*d4*" "*d5*" "*d6*" "*d7*" "*d8*"\n")
    
        write(f, "CELL_TYPES "*"1"*"\n")
        
        write(f, "12\n")          #四角のみ
    end     
    println("fin writing cells")
end

function write_result(out_file, outdir)
    temp1 = "CELL_DATA "*"1"*"\n"
    temp2 = "SCALARS "*"null"*" float\nLOOKUP_TABLE default\n"

    fff = outdir*"/"*out_file
    open(fff,"a") do f
        
        write(f, temp1)
        
        write(f, temp2)
        
        data = @sprintf("%7.7f", 1.0)
        write(f, data)
        write(f, "\n")
    
    end 
end

# main
main()

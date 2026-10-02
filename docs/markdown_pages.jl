"""
    front_matter_range(lines)

Return the range of `lines` holding the `---`-delimited YAML front matter that Weave
writes at the top of a markdown page, or an empty range if there is none.
"""
function front_matter_range(lines)
    (isempty(lines) || strip(first(lines)) != "---") && return 1:0
    closing = findnext(line -> strip(line) == "---", lines, 2)
    return closing === nothing ? (1:0) : (1:closing)
end

function unquote_yaml_scalar(value)
    value = strip(value)
    if length(value) >= 2 && first(value) == last(value) && first(value) in ('"', '\'')
        value = value[nextind(value, 1):prevind(value, lastindex(value))]
        value = replace(value, "\\\"" => "\"", "''" => "'")
    end
    return value
end

is_fence(line) = startswith(lstrip(line), "```") || startswith(lstrip(line), "~~~")

heading_level(line) = (m = match(r"^(#{1,6})(?:\s|$)", line)) === nothing ? 0 : length(m[1])

heading_text(line) = strip(rstrip(strip(line[(heading_level(line) + 1):end]), '#'))

"""
    benchmark_page(lines, fallback_title)

Turn the lines of a Weave-generated markdown file into `(title, body)` for Documenter.

The title is read from the `title:` key of the front matter wherever it appears (Weave
does not keep the keys in source order), falling back to `fallback_title`. The front
matter is dropped, as is a leading heading that just repeats the title. If the body uses
level-1 headings for its sections, every heading outside code blocks is demoted one level
so the page title stays the only level-1 heading on the page.
"""
function benchmark_page(lines, fallback_title)
    header = front_matter_range(lines)
    title = fallback_title
    for line in lines[header]
        key_value = split(line, ':'; limit = 2)
        if length(key_value) == 2 && strip(first(key_value)) == "title"
            title = unquote_yaml_scalar(last(key_value))
            break
        end
    end

    body = String.(lines[(last(header) + 1):end])
    first_content = findfirst(!isempty ∘ strip, body)
    if first_content !== nothing && heading_level(body[first_content]) > 0 &&
            lowercase(heading_text(body[first_content])) == lowercase(title)
        body = body[(first_content + 1):end]
    end

    in_code = false
    headings = Int[]
    for (i, line) in enumerate(body)
        if is_fence(line)
            in_code = !in_code
        elseif !in_code && heading_level(line) > 0
            push!(headings, i)
        end
    end
    if any(i -> heading_level(body[i]) == 1, headings)
        for i in headings
            heading_level(body[i]) < 6 && (body[i] = "#" * body[i])
        end
    end

    return title, body
end
